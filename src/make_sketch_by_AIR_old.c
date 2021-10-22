#include "kNN_buffer.h"
#include "Config.h"
#include "Sketch.h"
#include "ftr.h"

#ifndef FLIP_NUM
#define FLIP_NUM 1
#endif

/*
 ./make_Sketch database_reco0118.ftr test.sk 2 1
 argv[3]:(BP)    最適化方法   BP:1   QBP:2   PQBP:3
 とりあえず衝突法を用いて，QBPとPQBPに対してincrementalに最適化している．
 距離と順位付けの相関係数やAIRを用いる場合は，未対応．
 argv[3]:(AIR)   最適化方法   BP:1    <air + Exchenge>:2 <air + flip>:3 <air + swap>:4 
 argv[4]:   seed値
*/

int main(int argc, char *argv[])
{
	int i;
	// int kind = 1;
	// double db_score = 0;
//	char *database_file, *sample_database, *query_file, *answer_file, *sketch_file;
	char *sample_database, *query_file, *answer_file, *sketch_file;
	int seed, part, obj;

	if(argc != 9) {
		fprintf(stderr, "usage> %s sample_database query_file answer_file sketch_file lb partition 10*K seed\n", argv[0]);
		fprintf(stderr, "       lb:        Hamming -> 0, score_1 -> 1, score_inf -> 3\n");
		fprintf(stderr, "       partition: BP -> 1, QBP -> 2, PQBP -> 3\n");
		fprintf(stderr, "       10 * K:    precision for K = α%% < 5.0%% -> 10α, K for 90%% -> 90, K for 95%% -> 95\n");
		fprintf(stderr, "       seed:      any for random seed\n");
		return -1;
	}
//	fprintf(stderr, "K = %d, PJT_DIM=%d\n", K, PJT_DIM);
//	database_file = argv[1];
	sample_database = argv[1];
	query_file = 	argv[2];
	answer_file = 	argv[3];
	sketch_file =	argv[4];
	lb = 	   atoi(argv[5]);
	part =	   atoi(argv[6]);
	obj = 	   atoi(argv[7]);
	seed = 	   atoi(argv[8]);
	srandom(seed);

	if(!read_database(sample_database)) 	return -1;
//	fprintf(stderr, "read_database OK\n");
	if(!read_query(query_file)) 			return -1;
//	fprintf(stderr, "read_query OK\n");
	if(!read_answer(answer_file)) 			return -1;
//	fprintf(stderr, "read_answer OK\n");

	initialize(&pivot);
	pivot.partition = part;
	fprintf(stderr, "pivot.partition = %d (BP -> 1, QBP -> 2, PQBP -> 3)\n", pivot.partition);	
	fprintf(stderr, "objective function = %d (collision -> 0, precision for K = α%% -> α, K for 90%% -> 90, K for 95%% -> 95)\n", obj);	
	if(FTR_DIM != db_header->data_dim) {
		printf("FTR_DIM error\n");
		return -1;
	}
	fprintf(stderr, "get_sample: MAX_SAMPLESIZE = %d", MAX_SAMPLESIZE);
	get_sample(); // ピボット評価用のサンプル
	fprintf(stderr, " OK\n");
	median();//中央値計算
	fprintf(stderr, "median OK\n");
	clock_t start,end;
	time_t s_t, e_t;
	
	init_ind();
	flag = 1;
	start = clock(); s_t = time(NULL);
	fprintf(stderr, "QBP starts: NUM_TRIAL1 = %d\n", NUM_TRIAL1);
	if(NUM_TRIAL1 == 0) {
		random_QBP();
	} else {
		select_pivot_BPorQBP();
	}
	if(obj > 0) { // サンプルデータベースに対する検索精度を直接目的関数に用いて最適化する
		if(pivot.partition == 3 || FLIP_NUM == 1) {
//			optimize_pivot_by_precision_AIR(obj, EVAL_MODE); // 座標分割法およびFLIP＝1
		} else {
			fprintf(stderr, "FRIP = %d\n", FLIP_NUM);
			optimize_pivot_by_precision_AIR_flip(obj, EVAL_MODE, FLIP_NUM);
		}
	}
	end = clock(); e_t = time(NULL);
//	if(!read_database(database_file)) 		return -1;
//	fprintf(stderr, "read_database for test OK\n");
//	make_sketch();
	out_pivot(sketch_file);
//	out_sketch(sketch_file);
	fprintf(stderr, "TIME, %.2f, %5ld, TRIAL1, %d, N_0, %d, TRIAL2, %d, NUM_LS, %d, T_POWER, %.2f\n", (double)(end-start)/CLOCKS_PER_SEC, e_t - s_t, NUM_TRIAL1, N0, NUM_TRIAL2, NUM_LS, (double)T_POWER);
	printf(         "TIME, %.2f, %5ld, TRIAL1, %d, N_0, %d, TRIAL2, %d, NUM_LS, %d, T_POWER, %.2f\n", (double)(end-start)/CLOCKS_PER_SEC, e_t - s_t, NUM_TRIAL1, N0, NUM_TRIAL2, NUM_LS, (double)T_POWER);

	for(i = 0; i < db_header->data_num; i++) {
		free(database[i]);
		//free(database_temp[i]);
	}
	
	for(i = 0; i < MAX_SAMPLESIZE; i++)
			free(sample[i]);


//	print_result(seed, db_score,sc,(double)(end-start)/CLOCKS_PER_SEC);
	
	free(database);
	//free(database_temp);
	free(db_header);
	free(sub_fb);
	free(sketch);
	free(sample_sketch);
	for(i = 0; i < PJT_DIM; i++) {
		free(pivot.piv0[i]);
#if SKETCH_PARTITION == 0
		free(pivot.piv1[i]);
#endif
	}
	
	return 0;
}