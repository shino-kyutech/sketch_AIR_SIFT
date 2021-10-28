// FTR形式のデータファイルを間引いて（ratio毎に1個に間引いて）小さいデータファイルを作るプログラム
// （使用例）shrink_FTR ratio "output.ftr" "input_1.ftr" "input_2.ftr" ... 
// 第1パラメータ：縮小率（2 -> 2 分の 1，4 -> 4 分の 1）
// 第2パラメータ：出力ファイル
// 第3パラメータ 第4パラメータ ... ：入力ファイル1 入力ファイル2 ...
// 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "ftr.h"

int main(int argc, char *argv[])
{
	int ratio = atoi(argv[1]);
	char *output_ftr_filename = argv[2];
	int num_input_ftr_files = argc - 3;
	char **input_ftr_filename = argv + 3;

	fprintf(stderr, "output_ftr_file = %s\n", output_ftr_filename);
	fprintf(stderr, "input_ftr_files = ");
	for(int i = 0; i < num_input_ftr_files; i++) {
		fprintf(stderr, " %s", input_ftr_filename[i]);
	}
	fprintf(stderr, "\n");

	struct_multi_ftr *mf;
	mf = open_multi_ftr(num_input_ftr_files, input_ftr_filename, BLOCK_SIZE);
	int num_data = mf->num_data;
	fprintf(stderr, "open input files OK: total number of data = %d\n", num_data);

	ftr_header_type output_header;
	output_header = mf->hd[0]; 			// 入力ファイルの先頭のヘッダをコピー
	output_header.num_data = num_data / ratio; 	// データ数を入力の総数に書き換える

	FILE *fp;
	if((fp = fopen(output_ftr_filename, "wb"))  == NULL) {
		fprintf(stderr, "cannot open ftr_file: file name = %s\n", output_ftr_filename);
		exit(1);
	}
	fprintf(stderr, "open ftr_file OK: file name = %s\n", output_ftr_filename);
	fwrite(&output_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す

    struct_ftr_id *ftr_id_p;
	int i, written = 0;
	for(i = 0; i < num_data; i += ratio) {
		ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, &i, 0, 1);
		if(fwrite(ftr_id_p, sizeof(struct_ftr_id), 1, fp) != 1) {
			fprintf(stderr, "fwrite error: i = %d\n", i);
			exit(1);
		}
		if(!(written % 1000)) fprintf(stderr, "*");
        written++;
	}

    if(written != num_data / ratio) {
        output_header.num_data = written;
        fwrite(&output_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す
        fseek(fp, 0, SEEK_SET);
    }

	fclose(fp);

	fprintf(stderr, "written = %d\n", written);

	return 0;
}
