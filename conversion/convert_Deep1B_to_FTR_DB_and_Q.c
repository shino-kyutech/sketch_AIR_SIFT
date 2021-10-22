// Deep1Bの特徴データを.fbin（float）から.ftr（signed char）に変換するプログラム
// 元のデータにはデータ番号（ID）が付いていないので，読み込んだデータ順番に
// パラメタで指定するオフセットを加えたものをIDとする．
// 提供されたfbinファイルが3本ある．
// base.1B.fbin: 10億個のデータ
// base.10M.fbin:	1000万個のデータ（これが1Bのサブセットかどうかは調べてみる）
// base.350M.fbin:	3億5千万個のデータ（これと1Bや10Mとの関係も調べてみる）
// query.public.10K.fbin:	1万個の質問データ（これと一致するデータがあるか調べてみる）
// 後の取り扱いのために，10M個のデータ単位に分割する．
// 使用法：
// コンパイル：gcc -O3 -DBUFFER_SIZE=1000 -fopenmp convert_Deep1B_to_FTR_DB_and_Q.c -o convert_Deep1B_to_FTR_DB_and_Q
//   BUFFER_SIZE: まとめて読書きするレコード数（ただし，出力するレコード数を超えないようにする）
//   実行時パラメータ: fbin_file ftr_file from num
//   fbin_file = input file of fbin format
//   ftr_file  = output file ftr format
//   from      = data number of the 1st input data
//   num       = number of data to convert (0 for all to end)
// 元のfbin形式：先頭に32bit整数でデータ数と次元数が格納されている．
// それに続いて，次元数分のfloatが格納されている．
// 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

typedef signed char ftr_element_type; // 特徴データの要素の型
typedef ftr_element_type *ftr_type; // 特徴データベクトルの型
#define FTR_DIM 96 // 特徴データの次元数
#define DATA_NUM 1000000 // 出力するデータ総数（実際には，後で書換える）
#ifndef BUFFER_SIZE
#define BUFFER_SIZE 1000 // まとめて入力して，並列処理で変換し，まとめて出力する．
#endif
#ifndef CONVERT_FACTOR
#define CONVERT_FACTOR 255
#endif
#define CONVERT(v) ((int)((v) * CONVERT_FACTOR)) // float -> int の変換

// auxi_type（特徴データの付随情報）の型
typedef unsigned int auxi_type; // Deep1BデータではID（データ通し番号？）1B = 10^9 = 1G

// 特徴データファイルのヘッダ情報
typedef struct {
	int element_size;	// 各次元のサイズ（おそらく1(char)，または2(short)）
	int data_dim;		// 特徴次元数（image:64, music:96, colors:112, DeCAF:4096, Deep1B）
	int data_num;		// データ数（縮小するときは，後で変更する）
	int auxi_size;		// auxi_type（付随情報）のサイズ（Deep1Bでは，データ数が10憶，intに収まるので 4）
} ftr_header_type;

// 特徴データと付随情報（データID）をまとめた構造体
typedef struct {
	ftr_element_type ftr[FTR_DIM];
	auxi_type data_id;
} struct_ftr_id;

ftr_header_type db_header = {sizeof(ftr_element_type), FTR_DIM, DATA_NUM, sizeof(auxi_type)};

FILE *open_fbin(char *filename, int *num_data, int *dim);
void seek_fbin(FILE *f, int pos);	
int convert_vectors(struct_ftr_id ftr_id_buff[], float fbin_buff[][FTR_DIM], auxi_type from, int num, long *overflow);
int make_ftr_db(FILE *fbin, char *filename, int from, int num, int offset);

int main(int argc, char *argv[])
{
	int n, m;

	if(argc != 6) {
		fprintf(stderr, "usage> %s fbin_file ftr_file from num offset\n", argv[0]);
		fprintf(stderr, "fbin_file = input file of .fbin format\n");
		fprintf(stderr, "ftr_file  = output file .ftr format\n");
		fprintf(stderr, "from      = data number of the 1st input data\n");
		fprintf(stderr, "num       = number of data to convert (0 for all to end)\n");
		fprintf(stderr, "offset    = data id offset to add\n");
		return -1;
	}
	char *fbin_file = argv[1];
	char *ftr_file = argv[2];
	int from = atoi(argv[3]);
	int num = atoi(argv[4]);
	int offset = atoi(argv[5]);
	int num_data, dim;
	FILE *fbin = open_fbin(fbin_file, &num_data, &dim);
	if(dim != FTR_DIM) {
		fprintf(stderr, "dimension confliction, expected = %d, file = %d\n", FTR_DIM, dim);
		return -1;
	}
	if(num == 0) {
		num = num_data - from;
	} else if(num_data < from + num) {
		fprintf(stderr, "to small number of data in fbin file = %d, expected numbe of data = %d\n", num_data, from + num);
		return -1;
	}
	fprintf(stderr, "number of data in %s = %d, dimension of data = %d\n", fbin_file, num_data, dim);

	int num_db;
	if((num_db = make_ftr_db(fbin, ftr_file, from, num, offset)) != num) {
		fprintf(stderr, "cannot make %d records: number of records in output = %d\n", num, num_db);
		return -2;
	}
	
	fprintf(stderr, "total number of output data = %d\n", num_db);
	fclose(fbin);

	return 0;
}

FILE *open_fbin(char *filename, int *num_data, int *dim)
{
	FILE *f = fopen(filename, "rb");
	if(f == NULL) {
		fprintf(stderr, "cannot read open: file = %s\n", filename);
		exit(1);
	}
	if(fread(num_data, sizeof(int), 1, f) != 1) {
		fprintf(stderr, "read error (number of data): file = %s\n", filename);
		fclose(f);
		exit(1);
	}
	if(fread(dim, sizeof(int), 1, f) != 1) {
		fprintf(stderr, "read error (dim): file = %s\n", filename);
		fclose(f);
		exit(1);
	}
	return f;
}

void seek_fbin(FILE *f, int pos)
{
	long offset = pos * sizeof(float) * FTR_DIM + sizeof(int) * 2;
	if(fseek(f, offset, SEEK_SET) != 0) {
		fprintf(stderr, "fseek error\n");
		exit(1);
	}
}

// floatのvectorをsigned charのftrとdata_idに変換する 	
int convert_vectors(struct_ftr_id ftr_id_buff[], float fbin_buff[][FTR_DIM], auxi_type data_id, int num, long *overflow)
{
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(int r = 0; r < num; r++) {
		for(int i = 0; i < FTR_DIM; i++) {
			int val = CONVERT(fbin_buff[r][i]);
			if(val > 127) {
				val = 127;
				(*overflow)++;
			} else if(val < -127) {
				val = -127;
				(*overflow)++;
			}
			ftr_id_buff[r].ftr[i] = val;
		}
		auxi_type id = data_id + r;
		ftr_id_buff[r].data_id = id;
	}
}

// データベースを ftr 形式で作成する．
int make_ftr_db(FILE *fbin, char *filename, int from, int num, int offset)
{
	FILE *fp;
	int record_size = db_header.element_size * db_header.data_dim;
	ftr_type ftr = (ftr_type)malloc(record_size);
	int buffer_size = BUFFER_SIZE;
	if(buffer_size > num) buffer_size = num;
	long overflow = 0;

	if((fp = fopen(filename, "wb"))  == NULL) {
		fprintf(stderr, "cannot open ftr_file: file name = %s\n", filename);
		exit(1);
	}
	fprintf(stderr, "open ftr_file OK: file name = %s\n", filename);
	db_header.data_num = num;
	fwrite(&db_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す

//	static float fbin_buff[buffer_size][FTR_DIM];
//	static struct_ftr_id ftr_id_buff[buffer_size]; 
	static float (*fbin_buff)[FTR_DIM] = NULL;
	static struct_ftr_id *ftr_id_buff = NULL;
	if(fbin_buff == NULL) {
		fbin_buff = malloc(sizeof(float) * FTR_DIM * buffer_size);
	}
	if(ftr_id_buff == NULL) {
		ftr_id_buff = malloc(sizeof(struct_ftr_id) * buffer_size);
	}
	seek_fbin(fbin, from);
	int k = 0, data_id = from + offset;
	while(k < num) {
		int read_in = fread(fbin_buff, sizeof(float) * FTR_DIM, buffer_size, fbin);
		if(read_in == 0) break;
		convert_vectors(ftr_id_buff, fbin_buff, data_id, read_in, &overflow);
		if(fwrite(ftr_id_buff, sizeof(struct_ftr_id) * read_in, 1, fp) != 1) {
			fprintf(stderr, "fwrite error (ftr_id): data id = %d, record number = %d\n", data_id, k);
			exit(1);
		}
		k += read_in;
		data_id += read_in;
		fprintf(stderr, "*");
	}
	db_header.data_num = k;
	fseek(fp, 0L, SEEK_SET);
	fwrite(&db_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す

	fprintf(stderr, "overflow = %ld\n", overflow);
	fclose(fp);
	return k;
}
