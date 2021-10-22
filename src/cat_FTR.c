// FTR形式のデータファイルを連結するプログラム
// （使用例）cat_FTR "output.ftr" "input_1.ftr", ... , "input_n.ftr"
// 第1パラメータ：出力ファイル
// 第2パラメータ以降：入力ファイル
// 入力ファイルを連結して，出力ファイルを作成する．
// 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "ftr.h"

int main(int argc, char *argv[])
{
	int num_input_ftr_files = argc - 2;
	char **input_ftr_filename = argv + 2;
	char *output_ftr_filename = argv[1];

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
	output_header.num_data = num_data; 	// データ数を入力の総数に書き換える

	FILE *fp;
	if((fp = fopen(output_ftr_filename, "wb"))  == NULL) {
		fprintf(stderr, "cannot open ftr_file: file name = %s\n", output_ftr_filename);
		exit(1);
	}
	fprintf(stderr, "open ftr_file OK: file name = %s\n", output_ftr_filename);
	fwrite(&output_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す

	struct_ftr_id *ftr_id_p;
	int i;
	for(i = 0; i < num_data; i++) {
		ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, &i, 0, 1);
		if(fwrite(ftr_id_p, sizeof(struct_ftr_id), 1, fp) != 1) {
			fprintf(stderr, "fwrite error: i = %d\n", i);
			exit(1);
		}
		if(!(i % 10000000)) fprintf(stderr, "*");
	}
	fclose(fp);

	fprintf(stderr, "written = %d\n", i);

	return 0;
}
