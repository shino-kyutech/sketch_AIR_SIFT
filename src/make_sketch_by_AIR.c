#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "parm.h"
#include "config.h"
#include "bit_op.h"
#include "sketch.h"
#include "quick.h"
#include "pivot_selection.h"

int main(int argc, char *argv[])
{
	#ifdef _OPENMP
	#if defined(NUM_THREADS) && NUM_THREADS > 1
	omp_set_num_threads(NUM_THREADS);
	#endif
	int nt = omp_get_max_threads();
	#else
	int nt = 1;
	#endif
	fprintf(stderr, "number of threads = %d\n", nt);
	char *sample_dataset 	= SAMPLE_DATASET;
	char *query_file 		= QUERY_FILE;
	char *answer_file 		= ANSWER_FILE;
	char *pivot_file		= PIVOT_FILE;
	fprintf(stderr, "SEED = %d\n", SEED);
	srandom(SEED);
	#if defined(PARTITION_TYPE_QBP)
	int partition_type = QBP;
	#elif defined(PARTITION_TYPE_PQBP)
	int partition_type = PQBP;
	#else
	#error "PARTITION_TYPE_QBP or PARTITION_TYPE_PQBP should be defined."
	#endif
	struct_work_for_pivot_selection *work = new_work_for_pivot_selection(partition_type, sample_dataset, query_file, answer_file);
	clock_t start,end;
	time_t s_t, e_t;
	start = clock(); s_t = time(NULL);
	ftr_sample sample = {SAMPLE_SIZE_QBP, (ftr_type *)malloc(SAMPLE_SIZE_QBP * sizeof(ftr_type))};
	get_sample(work->ds_sample, &sample);
	select_pivot_QBP(work->pivot, work->med, &sample, work->ds_sample, nt);
	fprintf(stderr, "AIR starts: NUM_TRIAL_AIR = %d, TRUNCATE_AT = %d%%, EVAL_MODE = %d, NUM_FLIPS =%d", NUM_TRIAL_AIR, TRUNCATE_AT, EVAL_MODE, NUM_FLIPS);
	#ifdef VAR_FLIP
	fprintf(stderr, "VAR_FLIP\n");
	#else
	fprintf(stderr, "FIX_FLIP\n");
	#endif
	optimize_pivot_by_precision_AIR_flip(work);
	end = clock(); e_t = time(NULL);
	write_pivot(pivot_file, work->pivot);
	fprintf(stderr, "TIME, %.2f, %5ld, NUM_TRIAL_QBP, %d\n", (double)(end-start)/CLOCKS_PER_SEC, e_t - s_t, NUM_TRIAL_QBP);
	
	return 0;
}
