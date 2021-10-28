// FTR形式のデータファイルのヘッダ情報（要素サイズ，次元数，データ数）を表示するプログラム
// （使用例）info_FTR "filename.ftr"
// 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "ftr.h"

int main(int argc, char *argv[])
{
	char *input_ftr_filename = argv[1];
    ftr_header_type header;
    file_handle fh = open_ftr_file(input_ftr_filename, &header);
    print_header(stderr, input_ftr_filename, &header);
    CLOSE(fh);

	return 0;
}
