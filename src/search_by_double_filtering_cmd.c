#include <stdio.h>
#include <string.h>
#include "config.h"
#include "ftr.h"
#include "kNN_search.h"
#include "double_sketch.h"
#include "quick.h"
#include "bit_op.h"
#include "e_time.h"
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/resource.h>

int main(int argc, char *argv[])
{
//	#if NUM_K > 0
	int num_ftr_files = argc - 1;
	char **dataset_ftr_filename = argv + 1;
//	#endif
	char *expanded_pivot_file = EXPANDED_PIVOT_FILE;
	char *expanded_bucket_filename = EXPANDED_BUCKET_FILE;
	char *narrow_pivot_file = NARROW_PIVOT_FILE;
	char *narrow_bucket_filename = NARROW_BUCKET_FILE;
	char *query_ftr_filename = QUERY_FILE;
	char *query_2nd_ftr_filename = QUERY_2ND_FILE;
	char *query_3rd_ftr_filename = QUERY_3RD_FILE;
	char *qr_file[] = {query_ftr_filename, query_2nd_ftr_filename, query_3rd_ftr_filename};
//	#if NUM_K > 0
	char *answer_csv_filename = ANSWER_FILE;
	char *answer_2nd_csv_filename = ANSWER_2ND_FILE;
	char *answer_3rd_csv_filename = ANSWER_3RD_FILE;
	char *an_file[] = {answer_csv_filename, answer_2nd_csv_filename, answer_3rd_csv_filename};
	char *result_filename = RESULT_FILE;
	char temp_str[1000];
	char *result_filename2 = temp_str;
//	#endif
	int num_data;

	printf("NARROW_PJT_DIM = %d, PJT_DIM = %d: ", NARROW_PJT_DIM, PJT_DIM);
/*
	use_system("VmSize");
	struct_dataset *ds_query = read_dataset_n(1, &query_ftr_filename);
	int num_queries = ds_query->num_data;
	printf("read query file OK. the number of queries = %d: ", num_queries);
	use_system("VmSize");
	query_type *qr = (query_type *)malloc(sizeof(query_type) * num_queries);
	for(int i = 0; i < num_queries; i++) {
		qr[i] = (query_type) { i, ds_query->ftr_id[i].ftr};
	}

//	#if NUM_K > 0
	answer_type *correct_answer = read_correct_answer(answer_csv_filename, num_queries);
	printf("read correct answer OK. ");
	use_system("VmSize");
//	#endif
*/
	narrow_pivot_type *narrow_pivot = new_narrow_pivot(QBP);
	read_narrow_pivot(narrow_pivot_file, narrow_pivot);
	printf("read narrow pivot OK. ");
	use_system("VmSize");

	expanded_pivot_type *expanded_pivot = new_expanded_pivot(PQBP);
	read_expanded_pivot(expanded_pivot_file, expanded_pivot);
	printf("read expanded pivot OK. ");
	use_system("VmSize");

	struct_narrow_bucket *narrow_bucket_ds = read_narrow_bucket(narrow_bucket_filename);
	printf("read narrow bucket OK, ");
	num_data = narrow_bucket_ds->num_data;
	printf("number of data = %d: ", num_data);
	use_system("VmSize");

	struct_expanded_bucket *expanded_bucket_ds = read_compact_expanded_bucket(expanded_bucket_filename);
	printf("read compact expanded bucket OK, number of data = %d: ", expanded_bucket_ds->num_data);
	use_system("VmSize");
	
	if(num_data != expanded_bucket_ds->num_data) {
		fprintf(stderr, "narrow bucket is not compatible with expanded bucket: number of data in expanded bucket = %d\n", expanded_bucket_ds->num_data);
		return -1;
	}

	struct timespec tp01, tp02;
	clock_gettime(CLOCK_REALTIME, &tp01);
	// sort expanded sketches in the narrow sketch order
	fprintf(stderr, "sort expanded sketches in narrow sketch order, ");
	sort_expanded_sketches_narrow_sketch_order(expanded_bucket_ds->sk, narrow_bucket_ds->idx, num_data);
	clock_gettime(CLOCK_REALTIME, &tp02);
	fprintf(stderr, "OK. %.4lf: ", e_time(&tp01, &tp02));
	use_system("VmSize");

	int num_top_k = NUM_K; // NUM_K = kNN 検索で求める解の個数（ただし，-1: scoring のみ，0: filtering のみ）

	dataset_handle dh;

	#if defined(_OPENMP) && defined(NUM_THREADS)
		dh.num_threads = NUM_THREADS;
	#else
		dh.num_threads = 1;
	#endif
	
	dh.sorted = 0; // FTR is NOT sorted. Arrangement is as is for double filtering 
	dh.ftr_on = SECONDARY_MEMORY;
	dh.mf = (struct_multi_ftr **)malloc(sizeof(struct_multi_ftr *) * dh.num_threads);
	for(int t = 0; t < dh.num_threads; t++) {
		dh.mf[t] = open_multi_ftr(num_ftr_files, dataset_ftr_filename, BLOCK_SIZE);
	}
	dh.ds = NULL;
	if(num_data != dh.mf[0]->num_data) {
		fprintf(stderr, "bucket (filename = %s, num_data = %d) is not compatible with datasets (filename = %s, ... , num_data = %d)\n", narrow_bucket_filename, num_data, dataset_ftr_filename[0], dh.mf[0]->num_data);
		return -1;
	}

	#ifdef _OPENMP
	#if defined(NUM_THREADS) && NUM_THREADS > 1
	omp_set_num_threads(NUM_THREADS);
	#endif
	#endif

	struct_query_narrow_sketch *query_narrow_sketch = (struct_query_narrow_sketch *)malloc(sizeof(struct_query_narrow_sketch));
	struct_query_expanded_sketch *query_expanded_sketch = (struct_query_expanded_sketch *)malloc(sizeof(struct_query_expanded_sketch));
	struct_que_c2_n *que = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n));
	printf("allocate que for sketch enumeration OK. ");
	use_system("VmSize");
//	struct_que *que = (struct_que *)malloc(sizeof(struct_que));
	
	dist_type *score = (dist_type *)malloc(sizeof(dist_type) * num_data);
	int *idx = (int *)malloc(sizeof(int) * num_data);
	#ifdef SELECT_K_BY_QUICK_SELECT_K_PARA_WORK
	int *idx_temp = (int *)malloc(sizeof(int) * num_data);
	#endif
	#ifdef _OPENMP
	#if defined(NUM_THREADS) && NUM_THREADS > 1
	#ifdef SELECT_K_BY_QUICK_SELECT_K_PARA_WORK
	work_select_para *wp = new_work_select_para(NUM_THREADS);
	#endif
	#endif
	#endif
	printf("allocate memories for filtering OK. ");
	use_system("VmSize");

	int *data_num = (int *)malloc(sizeof(int) * num_data);
	int *data_num_temp = (int *)malloc(sizeof(int) * num_data);

	int num_query_files = sizeof(qr_file) / sizeof(qr_file[0]);
	query_type *qr[num_query_files];
	answer_type *correct_answer[num_query_files];
	int num_queries;
	for(int m = 0; m < num_query_files; m++) {
		if(strcmp(qr_file[m], "NONE") == 0) {
			num_query_files = m;
			break;
		}
		use_system("VmSize");
		struct_dataset *ds_query = read_dataset_n(1, &qr_file[m]);
		num_queries = ds_query->num_data;
		printf("read query file (%s) OK. the number of queries = %d: ", qr_file[m], num_queries);
		use_system("VmSize");
		qr[m] = (query_type *)malloc(sizeof(query_type) * num_queries);
		for(int i = 0; i < num_queries; i++) {
			qr[m][i] = (query_type) { i, ds_query->ftr_id[i].ftr};
		}
		// free_dataset(ds_query);
		correct_answer[m] = read_correct_answer(an_file[m], num_queries);
		printf("read correct answer (%s) OK. ", an_file[m]);
		use_system("VmSize");
	}

	#ifdef NUM_Q
	num_queries = NUM_Q;
	#endif
	kNN_buffer **top_k = (kNN_buffer **)malloc(sizeof(kNN_buffer *) * num_queries);
	for(int i = 0; i < num_queries; i++) {
		top_k[i] = NULL;
	}

    char cmd[1000], *cmdline;
    while(1) {
        fprintf(stderr, "nc1 (ppm) ? ");
        cmdline = fgets(cmd, 900, stdin);
        if(cmdline == NULL) break;
        double nc1 = atof(cmd);
        fprintf(stderr, "nc1 = %lf\n", nc1);
        if(nc1 < 0) break;
        fprintf(stderr, "nc2 (ppm) ? ");
        cmdline = fgets(cmd, 900, stdin);
        if(cmdline == NULL) break;
        double nc2 = atof(cmd);
        fprintf(stderr, "nc2 = %lf\n", nc2);
        fprintf(stderr, "OK -> <enter>, reset -> -1");
        cmdline = fgets(cmd, 900, stdin);
        if(cmdline == NULL) break;
        if(cmd[0] != '\n') continue;

		int num_candidates_1st = nc1 * 0.000001 * num_data;
		int num_candidates_2nd = nc2 * 0.000001 * num_data;
        fprintf(stderr, "num_data = %d, NUM_CANDIDATES_1st = %.2lf (ppm), num_candidates_1st = %d, NUM_CANDIDATES_2nd = %.2lf (ppm), num_candidates_2nd = %d\n", num_data, nc1, num_candidates_1st, nc2, num_candidates_2nd);

        for(int m = 0; m < num_query_files; m++) {
            use_system("VmSize");

            struct timespec tp1, tp2, tp3, tp4, tp5, tp6, tp7;
            clock_gettime(CLOCK_REALTIME, &tp1);
            double e_time_1st = 0, e_time_score = 0, e_time_2nd = 0, e_time_kNN = 0;
                
            for(int i = 0; i < num_queries; i++) {
                clock_gettime(CLOCK_REALTIME, &tp3);
                // 1st filtering by sketch enumeration
                set_query_narrow_sketch_p(query_narrow_sketch, &qr[m][i], narrow_pivot, SCORE_P_1ST / 10.0);
                filtering_by_sketch_enumeration_c2_n(query_narrow_sketch, narrow_bucket_ds, que, data_num, num_candidates_1st);
                clock_gettime(CLOCK_REALTIME, &tp4);
                e_time_1st += e_time(&tp3, &tp4);

                // 2nd filtering by sequential filtering using expanded sketches
                set_query_expanded_sketch_p(query_expanded_sketch, &qr[m][i], expanded_pivot, SCORE_P_2ND / 10.0);
                #if defined(_OPENMP) && defined(NUM_THREADS) && NUM_THREADS > 1 
                omp_set_num_threads(NUM_THREADS);
                #pragma omp parallel for
                #endif
                for(int j = 0; j < num_candidates_1st; j++) {
                    #ifdef SELECT_K_BY_QUICK_SELECT_K_PARA_WORK
                    idx_temp[j] = j;
                    #else
                    idx[j] = j;
                    #endif
                    score[j] = expanded_priority(expanded_bucket_ds->sk[data_num[j]], query_expanded_sketch); // sorting
                }
                clock_gettime(CLOCK_REALTIME, &tp5);
                e_time_score += e_time(&tp4, &tp5);

                // top-k selection for 2nd filtering
                int k;
                static int fst = 1;
                #if defined(SELECT_K_BY_QUICK_SELECT_K)
                if(fst) { fprintf(stderr, "quick_select_k\n"); fst = 0; }
                quick_select_k(idx, score, 0, num_candidates_1st - 1, num_candidates_2nd);
                for(int k = 0; k < num_candidates_2nd; k++) {
                    data_num_temp[k] = narrow_bucket_ds->idx[data_num[idx[k]]]; // sorting
                }
                #elif defined(SELECT_K_BY_QUICK_SELECT_K_PARA_WORK)
                if(fst) { fprintf(stderr, "quick_select_k_para_work\n"); fst = 0; }
                k = quick_select_k_para_work(idx_temp, score, num_candidates_1st, num_candidates_2nd, wp);
                int k2 = get_top_k_from_work_para(idx, idx_temp, wp);
                if(k != k2) fprintf(stderr, "k = %d, k2 = %d\n", k, k2);
                for(k = 0; k < k2; k++) {
                    data_num_temp[k] = narrow_bucket_ds->idx[data_num[idx[k]]]; // sorting
                }
                #else
                if(fst) { fprintf(stderr, "quick_sort\n"); fst = 0; }
                quick_sort(idx, score, 0, num_candidates_1st - 1);
                for(int k = 0; k < num_candidates_2nd; k++) {
                    data_num_temp[k] = narrow_bucket_ds->idx[data_num[idx[k]]]; // sorting
                }
                #endif
                k = num_candidates_2nd;
                clock_gettime(CLOCK_REALTIME, &tp6);
                e_time_2nd += e_time(&tp4, &tp6);

                // k-NN search from candidates
                if(top_k[i] == NULL) {
                    top_k[i] = new_kNN_buffer(num_top_k);
                }
                if(num_top_k == 1) {
                    search_NN(&dh, &qr[m][i], k, data_num_temp, top_k[i]);
                } else {
                    search_kNN(&dh, &qr[m][i], k, data_num_temp, top_k[i]);
                }
                
                clock_gettime(CLOCK_REALTIME, &tp7);
                e_time_kNN += e_time(&tp6, &tp7);

                #ifdef PRINTOUT_EACH_QUERY
                #if NUM_Q <= 100
                fprintf(stderr, "%.4lf, %.4lf, %d, %d, %d\n", e_time(&tp6, &tp7), e_time(&tp3, &tp7), i, top_k[i]->buff[0].data_num, top_k[i]->buff[0].dist);
                #endif
                #endif
            }

            clock_gettime(CLOCK_REALTIME,&tp2);

            out_result(result_filename, num_queries, correct_answer[m], top_k);
            strcpy(result_filename2, result_filename);
            result_filename2 = strtok(result_filename2, ".");
            result_filename2 = strcat(result_filename2, "_result.csv");
            double recall = out_result_double(result_filename2, NARROW_PJT_DIM, EXPANDED_PJT_DIM, FTR_DIM, e_time_1st, e_time_score, e_time_2nd, e_time_kNN, e_time(&tp1, &tp2), nc1, nc2,
                                                num_queries, correct_answer[m], top_k, qr_file[m]);
            printf("1st, %.4lf, score, %.4lf, 2nd, %.4lf, kNN, %.4lf, total, %.4lf, recall, %.1lf\n", e_time_1st, e_time_score, e_time_2nd, e_time_kNN, e_time(&tp1, &tp2), recall);
		}
	}

	return 0;
}
