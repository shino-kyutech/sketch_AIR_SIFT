#include "parm.h"
#include "bit_op.h"
#include "kNN_search.h"
#include "pivot_selection.h"
#include "sketch.h"
#include "quick.h"
#include <omp.h>

// dataset ds からsample->num_samples個のサンプルを取得して，配列sample->ftrに格納する
// ただし，sample が未割当（NULL）のときは，動的に割り当てる
void get_sample(struct_dataset *ds, ftr_sample *sample)
{
	int sample_number[sample->num_samples];
	if(sample->ftr == NULL) {
		sample->ftr = (ftr_type *)malloc(sizeof(ftr_type) * sample->num_samples);
	}
	for(int i = 0; i < sample->num_samples;){
		int x = random() % ds->num_data;
		int j;
		for(j = 0; j < i; j++) {
			if(x == sample_number[j]) break; // 同じデータは使わないようにする
		}
		if(j < i) continue;
		sample->ftr[i] = ds->ftr_id[x].ftr;
		sample_number[i++] = x;
	}
}

// データセットの中央値配列 med を求める
// 配列 med を動的に確保し，各次元 j = 0, ... , FTR_DIM - 1 に対する中央値 med[j] を求める
// ただし，すべてのデータを調べると時間がかかるので，500 個毎に一つのサンプルで求める
// DeCAF のときは，非負 8-bit 整数なので，0 ～ 255 の範囲の中央値を求める．
// そのときは，DECAF を #define すること．
// そうでないときは，符号付 8-bit 整数なので，-128 から 127 の範囲の中央値を求める．
ftr_type get_median(struct_dataset *ds)
{
	long cnt[256];
	int i, j, s, t;
	ftr_type med = (ftr_type)malloc(FTR_DIM * sizeof(ftr_element_type));

	for(j = 0; j < FTR_DIM; j++) {
		for(i = 0; i < 256; i++) {
			cnt[i] = 0;
		}
		for(i = t = 0; i < ds->num_data; i += 500) {
			#if defined(DECAF)
			cnt[ds->ftr_id[i].ftr[j]]++; t++;
			#else
			cnt[ds->ftr_id[i].ftr[j] + 127]++; t++;
			#endif
		}
		for(i = s = 0; (i < 256) && (s < t / 2); i++) {
			s += cnt[i];
		}
		#if defined(DECAF)
		med[j] = i;
		#else
		med[j] = i - 127;
		#endif
	}
	return med;
}

#ifndef NUM_TRIAL_QBP
#define NUM_TRIAL_QBP 1
#endif

// 乱択したデータを2値量子化した点を用いてピボットを作る．
// とりあえずのバージョン
// 最小衝突法は，素朴にサンプルの衝突を計算する版を用いる．衝突計算の高速化版は未実装（後回し）．
// pivot_type *pivot;
// ftr_type median;
// ftr_type sample[];	// sample data points
// struct_dataset *ds;	// dataset
// int nt;				// Number of threads
// Macros:
//    	SAMPLE_SIZE_QBP			= 衝突計算に用いるサンプル数
// 		FTR_DIM, PJT_DIM		= 特徴データの次元数，スケッチ幅（w）
// 		NUM_TRIAL_QBP			= 各次元でピボット選択をする試行回数（衝突が最小のものを選択する）
// 		SAMPLE_SIZE_RADIUS		= ピボットの半径を計算するときのサンプル数（ピボットの中心点とサンプル間の距離を中央値を求める）
void select_pivot_QBP(pivot_type *pivot, ftr_type median, ftr_sample *sample, struct_dataset *ds, int nt)
{
	static dist_type *sample_dist = NULL;
	static int *idx = NULL;
	static sketch_type *sample_sketch = NULL;
	if(sample_dist == NULL)
		sample_dist = (dist_type *)malloc(sizeof(dist_type) * SAMPLE_SIZE_QBP);
	if(idx == NULL)
		idx = (int *)malloc(sizeof(int) * SAMPLE_SIZE_QBP);
	if(sample_sketch == NULL)
		sample_sketch = (sketch_type *)malloc(sizeof(sketch_type) * SAMPLE_SIZE_QBP);
	dist_type min_r = 0;				// それまでに求めた最良のピボットの半径
	ftr_element_type min_p[FTR_DIM];	// 								 のセンター点

	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(int i = 0; i < SAMPLE_SIZE_QBP; i++) {
		#ifndef EXPANDED_SKETCH
		sample_sketch[i] = 0;
		#else
		for(int j = 0; j < SKETCH_SIZE; j++) {
			sample_sketch[i][j] = 0;
		}
		#endif
	}
	
	// QBP で使用する乱数がget_threshold_kやquick_select内部で使用する乱数と干渉するのを避けるために，ここでまとめて乱数を取得しておく
	int rdm[NUM_TRIAL_QBP * PJT_DIM];
	for(int i = 0; i < NUM_TRIAL_QBP * PJT_DIM; i++) rdm[i] = random() % ds->num_data;
	int ir = 0;

	int dim;
	double scmin;
	for(dim = 0; dim < PJT_DIM; dim++) {
		scmin = DBL_MAX;
		for(int t = 0; t < NUM_TRIAL_QBP; t++){	// 試行回数
			int c = rdm[ir++];
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < FTR_DIM; i++) {
				#ifndef PARTITION_TYPE_PQBP
				pivot->p[dim][i] = ds->ftr_id[c].ftr[i] < median[i] ? FTR_MIN : FTR_MAX;
				#else
				if(i < PART_START(dim) || i >= PART_START(dim) + PART_DIM(dim)) {
					pivot->p[dim][i] = 0; // 使用しないところはすべて0にしておく．
				} else {
					pivot->p[dim][i] = ds->ftr_id[c].ftr[i] < median[i] ? FTR_MIN : FTR_MAX;
				}
			#endif
			}
			int ss = SAMPLE_SIZE_RADIUS;
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < ss; i++) {
				int j = random() % ss;
				#ifndef PARTITION_TYPE_PQBP
				sample_dist[i] = DISTANCE(pivot->p[dim], sample->ftr[j], FTR_DIM);
				#else
				sample_dist[i] = PART_DISTANCE(pivot->p[dim], sample->ftr[j], PART_START(dim), PART_DIM(dim));
				#endif
				idx[i] = i;
			}
			if(ss > 100) {
				pivot->r[dim] = get_threshold_k(idx, sample_dist, ss, ss / 2, nt);  // 半径計算(get_thresholdを用いて中央値を求めている)
//				pivot->r[dim] = get_threshold_k(idx, sample_dist, ss, 3 * ss / 9 + random() % (2 * ss / 9), nt);  // 半径計算(中央値より少しずらしたものを半径にする)
			} else {
				insertion_sort(sample_dist, ss); // 半径計算に使用するサンプル数が少ないときは insertion_sort を用いて中央値を計算する
				pivot->r[dim] = sample_dist[ss / 2];
			}
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < sample->num_samples; i++) {
				#ifndef EXPANDED_SKETCH
				data_to_sketch_1bit(sample->ftr[i], pivot, dim, &sample_sketch[i]);
				#else
				data_to_sketch_1bit(sample->ftr[i], pivot, dim, sample_sketch[i]);
//				if(dim == 3 && i < 10) {
//					printf("i = %2d:", i); print_bin_expanded(sample_sketch[i], SKETCH_SIZE); printf("\n"); getchar();
//				}
				#endif
			}
			double sc = collision(SAMPLE_SIZE_QBP, sample_sketch);
			if(sc < scmin) {
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(int i = 0; i < FTR_DIM; i++) {
					min_p[i] = pivot->p[dim][i];
				}
				scmin = sc;
				min_r = pivot->r[dim];
			}
		}

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i = 0; i < FTR_DIM; i++) {
			pivot->p[dim][i] = min_p[i]; //best pivot
		}
		pivot->r[dim] = min_r;           //best radian
		fprintf(stderr,"DIM%3d fin ... = %12.0f\n", dim, scmin);
		if(scmin <= 1.0) {
			fprintf(stderr, "replace sample\n");
			get_sample(ds, sample);
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < SAMPLE_SIZE_QBP; i++) {  // 新しいsampleのsketchを作り直す
				for(int j = 0; j <= dim; j++) {
					#ifndef EXPANDED_SKETCH
					data_to_sketch_1bit(sample->ftr[i], pivot, j, &sample_sketch[i]);
					#else
					data_to_sketch_1bit(sample->ftr[i], pivot, j, sample_sketch[i]);
					#endif
				}
			}
		} else {
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < SAMPLE_SIZE_QBP; i++) { // Bestなpivotでsketchを作り直す（Trial中にsketchを書き換えているため）
				#ifndef EXPANDED_SKETCH
				data_to_sketch_1bit(sample->ftr[i], pivot, dim, &sample_sketch[i]);
				#else
				data_to_sketch_1bit(sample->ftr[i], pivot, dim, sample_sketch[i]);
				#endif
			}
		}
	}
	fprintf(stderr,"DIM%3d fin ... = %12.0f, %12.4e\n", dim, scmin, scmin / ((double)SAMPLE_SIZE_QBP * (SAMPLE_SIZE_QBP - 1) / 2));

	fprintf(stderr, "XXX\n");
	// サンプルデータの先頭分の100個のスケッチを作って，csv形式で書き出す（stdout）（EXPANDEDは未対応）
/*
	#ifndef EXPANDED_SKETCH
	sketch_type sk;
	for(int i = 0; i < 100; i++) {
		printf("i = %d\n", i);
		for(int j = 0; j < PJT_DIM; j++) {
			data_to_sketch_1bit(ds->ftr_id[i].ftr, pivot, j, &sk);
		}
		for(int j = 0; j < FTR_DIM; j++) {
			printf("%4d,", ds->ftr_id[i].ftr[j]);
		}
		printf("\n%lu,'", (long)sk);
		print_bin_long((long)sk);
		printf("\n");
	}
	#endif
*/
}

// とりあえず，乱択したデータを2値量子化した点を用いてピボットを作る．
// 最小衝突法は，用いない．
// 半径を求めるときにサンプル数を少なくして，毎回サンプルを選び直して，中央値を求める．
// サンプルが少ないと，中央値はばらつくので，そのばらつきを利用して，アンバランスな分割も混ざるようにする．
void select_pivot_random_QBP(pivot_type *pivot, ftr_type median, struct_dataset *ds)
{
	static dist_type sample_dist[SAMPLE_SIZE_RADIUS];

	// QBP で使用する乱数がget_threshold_kやquick_select内部で使用する乱数と干渉するのを避けるために，ここでまとめて乱数を取得しておく
	int rdm[PJT_DIM];
	for(int i = 0; i < PJT_DIM; i++) rdm[i] = random() % ds->num_data;

	int dim;
	for(dim = 0; dim < PJT_DIM; dim++) {
		int c = rdm[dim];
		for(int i = 0; i < FTR_DIM; i++) {
			#ifndef PARTITION_TYPE_PQBP
			pivot->p[dim][i] = ds->ftr_id[c].ftr[i] < median[i] ? FTR_MIN : FTR_MAX;
			#else
			if(i < PART_START(dim) || i >= PART_START(dim) + PART_DIM(dim)) {
				pivot->p[dim][i] = 0; // 使用しないところはすべて0にしておく．
			} else {
				pivot->p[dim][i] = ds->ftr_id[c].ftr[i] < median[i] ? FTR_MIN : FTR_MAX;
			}
			#endif
		}
		// 半径計算
		for(int i = 0; i < SAMPLE_SIZE_RADIUS; i++) {
			int x = random() % ds->num_data;
			#ifndef PARTITION_TYPE_PQBP
			sample_dist[i] = DISTANCE(pivot->p[dim], ds->ftr_id[x].ftr, FTR_DIM);
			#else
			sample_dist[i] = PART_DISTANCE(pivot->p[dim], ds->ftr_id[x].ftr, PART_START(dim), PART_DIM(dim));
			#endif
		}
		insertion_sort(sample_dist, SAMPLE_SIZE_RADIUS);
		pivot->r[dim] = sample_dist[random() % SAMPLE_SIZE_RADIUS]; 
	}
}

#ifndef EXPANDED_SKETCH

static int comp_sketch_2(const void *a, const void *b) {

	if(*(sketch_type *)a < *(sketch_type *)b)
		return -1;
	else if(*(sketch_type *)a == *(sketch_type *)b)
		return 0;
	else
		return 1;
	
}

#else

static int comp_sketch_2(const void *a, const void *b) 
{
	sketch_type *x = (sketch_type *)a, *y = (sketch_type *)b;
	for(int j = 0; j < SKETCH_SIZE; j++) {
		if((*x)[j] < (*y)[j])
			return -1;
		else if((*x)[j] == (*y)[j])
			continue;
		else
			return 1;
	}
	return 0;
}

#endif

// sk[0], ... , sk[SAMPLE_SIZE - 1] 内の衝突の組合せ数を求める．
// 値がオーバーフローしないように，関数の返り値は double にしている．
// sk の要素数は，マクロのSAMPLE_SIZE．ただし，関数の引数に変更した方がよいかも．
// 現状の衝突の求め方は，素朴なものなので，高速化版の実装を復元した方がよい．
double collision(int num_sk, sketch_type sk[])
{
	double n;       // 衝突している組の個数
 	double score = 0;       // スコア
	static sketch_type *temp = NULL;
	if(temp == NULL) {
		temp = (sketch_type *)malloc(sizeof(sketch_type) * num_sk);
	}
	for(int i = 0; i < num_sk; i++) {
		#ifndef EXPANDED_SKETCH
		temp[i] = sk[i];
		#else
		for(int j = 0; j < SKETCH_SIZE; j++) {
			temp[i][j] = sk[i][j];
		}
		#endif
	}
	qsort(temp, num_sk, sizeof(sketch_type), comp_sketch_2); // skのコピーをソート
	for (int i = 0; i < num_sk - 1; i++) {
		n = 0;
		#ifdef EXPANDED_SKETCH
		while(i < num_sk - 1 && comp_sketch_2(&temp[i], &temp[i + 1]) == 0) {
			n++;
			i++;
		}
		#else
		while(i < num_sk - 1 && temp[i] == temp[i + 1]) {
			n++;
			i++;
		}
		#endif
		if(n > 0) {
			n++;
			score += n * (n - 1) / 2;
		}
	}
	return score;
}

// wide スケッチ用は，まずは，QBP のみで 
// optimize_pivot_by_precision // LS and AIR

#ifdef USE_AIR
// create work space for pivot selection by LS(local search) or AIR

// 目的関数の指定（精度）TRUNCATE_AT (1 ～ 89) (パーミル ‰) = K で候補数を打ち切って検索したときの精度
// ただし，w (PJT_DIM) >= 64 のときは，精度が容易に高くなるので，ppm で与える
//                    （候補数）90 <= obj < 100 → 精度が obj (%) になる候補数 K(%) (未調整)
// ただし，現状では上の精度でのみデバッグ・調整中
// 評価の集計方法（EVAL_MODE: 0 -> 正解数の割合，1 -> 正解の重み付きスコア合計）
// NUM_FLIPS:  1回の試行で同時にフリップ（0 と 255 を反転）する次元数
void optimize_pivot_by_precision_AIR_flip(struct_work_for_pivot_selection *work)
{
	// 【注意】とりあえず，座標分割法には未対応なので，SKETCH_PARTITION = 3 では使用しないこと（ただいま対応のための変更中）
	// 質問集合に対して，再サンプリングの最終サイズを半分で止めるようにする．#define USE_HALF 
	// USE_HALF は精度を目的関数とするときのみ，つまり，1 <= obj <= 89 にのみ対応

    pivot_type *pivot = work->pivot;
    struct_dataset *ds_sample = work->ds_sample;
    struct_dataset *ds_query = work->ds_query;
// int obj, int mode, int num_flips,
	int n0 = N0;
	int reuse = REUSE;
	int resampling_size = 10;
	double Tr;
	int num_flips = NUM_FLIPS;
	#ifdef VAR_FLIP
		#ifndef REPLACE_WHOLE
			int f_max = NUM_FLIPS; // フリップの最大数
		#endif
	#endif
	int not_improved = 0, final_check = 0;
	#ifndef CONVERGENCE_CHECK
		#define CONVERGENCE_CHECK 0
	#endif
	int num_queries = ds_query->num_data;
	int max_num_queries = num_queries;
	int num_accept = 0;
	#ifdef USE_HALF
		max_num_queries /= 2;
	#endif
	#ifdef REPLACE_WHOLE
		#ifndef REPLACE_TRIAL
			#define REPLACE_TRIAL 1
		#endif
		int j;
		int candidate[REPLACE_TRIAL];
		static ftr_type piv_candidate[REPLACE_TRIAL] = {NULL};
		static ftr_type piv_temp = NULL;
		int r_count;
		int piv_diff[REPLACE_TRIAL], piv_diff_min, piv_diff_min_n = 0;
	#endif
	ftr_sample sample_for_radius = {SAMPLE_SIZE_RADIUS, (ftr_type *)malloc(SAMPLE_SIZE_RADIUS * sizeof(ftr_type))};
	get_sample(ds_sample, &sample_for_radius);
	dist_type rad_temp;
	int dim, axi[FTR_DIM];
	int /* multi_K[1000], */ accept;
	double eval, eval_best, prec, prec_best;

//	fprintf(stderr, "make_sketch ... ");
	make_sketch(work);
//	fprintf(stderr, "OK\nmake_idx_bkt ... ");
	make_idx_bkt(work->sketch, work->ds_sample->num_data, &work->bucket->idx, &work->bucket->bkt);
//	fprintf(stderr, "OK\n");
	for(int i = 0; i < FTR_DIM; i++) {
		axi[i] = i;
	}
	int obj = TRUNCATE_AT, mode = EVAL_MODE;
	make_query_sketch_for_resample(max_num_queries, work);
	eval_best = EVAL_PRECISION(obj, max_num_queries, mode, &prec_best, work);
	#if PJT_DIM < 64
	printf("(1) Precision for K = %.1lf%% = %10.8lf%%, prec = %10.8lf, reuse = %d\n", (double)obj / 10, eval_best, prec_best, reuse);
	#else
	printf("(1) Precision for K = %d ppm = %10.8lf%%, prec = %10.8lf, reuse = %d\n", obj, eval_best, prec_best, reuse);
	#endif

	// なぜか，ここでresample_queryを1000回行っていたが，理由が不明（たぶん不要なのでコメントアウト）
	// for(i = 0; i < 1000; i++) {
	// 	resample_query(max_num_queries, work);
	// }
	
	for(int i = 0, reuse_c = 0; i < NUM_TRIAL_AIR + NUM_TRIAL_LS + CONVERGENCE_CHECK; i++, reuse_c = (reuse_c + 1) < reuse ? reuse_c + 1 : 0) {
		if(i == NUM_TRIAL_AIR) { // DO Local Search using whole queries
			resampling_size = max_num_queries;
			fprintf(stderr, "Goto local search: trials = %d\n", NUM_TRIAL_LS);
			printf("Goto local search: trials = %d\n", NUM_TRIAL_LS);
			make_query_sketch_for_resample(resampling_size, work);
			// eval_best = precision_resample_cursor_2_n(obj, resampling_size, mode, &prec, work);
			eval_best = EVAL_PRECISION(obj, resampling_size, mode, &prec_best, work);
			#if PJT_DIM < 64
			printf("(2) Precision for K = %.1lf%% = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", (double)obj / 10, eval_best, prec_best, resampling_size, i);
			#else
			printf("(2) Precision for K = %d ppm = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", obj, eval_best, prec_best, resampling_size, i);
			#endif
		}
		if(i == NUM_TRIAL_AIR + NUM_TRIAL_LS) {
			fprintf(stderr, "Goto convergence check: trials = %d\n", CONVERGENCE_CHECK);
			printf("Goto convergence check: trials = %d\n", CONVERGENCE_CHECK);
			final_check = 1;
		} else if(i < NUM_TRIAL_AIR && reuse_c == 0) {
			// AIR では，REUSE回毎に再サンプリングを行う．
			int Ti = i + REUSE;
			if(Ti > NUM_TRIAL_AIR) Ti = NUM_TRIAL_AIR - 1;
			Tr = 1 - (double)Ti / (NUM_TRIAL_AIR == 0 ? 1 : NUM_TRIAL_AIR); // Tr = remain; 1 - Tr = progress;
			if(T_POWER >= 1.0) {
				Tr = pow(Tr, T_POWER);
			} else {
				Tr = Tr * pow(T_POWER, 1 - Tr); // remain^aplha * beta^progress; alpha = 1, beta = T_POWER;
			}
			resampling_size = (double)max_num_queries / ((double)(max_num_queries - n0) / n0 * Tr * Tr + 1);
			reuse_c = 0;

//			printf("Resampling: size = %d\n", resampling_size);
			resample_query(resampling_size, work);
//			for(int d = 0; d < resampling_size; d++) {
//				printf("query_num[%1d] = %d\n", d, work->query_sketch[d].query.query_num);
//			}
//			getchar();
			make_query_sketch_for_resample(resampling_size, work);
			// eval_best = precision_resample_cursor_2_n(obj, resampling_size, mode, &prec, work);
			eval_best = EVAL_PRECISION(obj, resampling_size, mode, &prec_best, work);
			#if PJT_DIM < 64
			printf("(3) Precision for K = %.1lf%% = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", (double)obj / 10, eval_best, prec_best, resampling_size, i);
			#else
			printf("(3) Precision for K = %d ppm = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", obj, eval_best, prec_best, resampling_size, i);
			#endif
		} else {
			// LS（CONVERGENCE_CHECK）では，再サンプリングは行わない．
		}
		// ピボット（dim = 0, ... , PJT_DIM)と次元 (axi = 0, ... , FTR_DIM) を乱択
		// pivot.piv0[dim][axi] の値を反転（フリップ）（0 <==> 255 (DeCAF)) or (-127 <==> 127 (Deep1B))
		// 半径を更新
		if(final_check) {
			if(final_check == 1) {
				dim = 0;
				axi[0] = 0;
				num_flips = 1;
				final_check = 2;
			} else if(final_check == 2) {
				axi[0]++;
				if(axi[0] >= FTR_DIM) {
					axi[0] = 0;
					dim++;
					if(dim >= PJT_DIM) {
						dim = 0;
					}
				}
			}
			fprintf(stderr, "Convergence check (i = %7d): dim = %2d, axi = %2d, not_improved = %4d\r", i - NUM_TRIAL_AIR - NUM_TRIAL_LS, dim, axi[0], not_improved);
		} else {
			// fprintf(stderr, "Try FLIP\n");
			dim = random() % PJT_DIM;
			#ifdef REPLACE_WHOLE
				// ピボットを一つ(dim次元のもの)を全部取り換える．REPLACE_TRIAL個の候補から現在のピボットとの違いが最小のものを選ぶ．
				if(piv_candidate[0] == NULL) {
					for(j = 0; j < REPLACE_TRIAL; j++) { 
						piv_candidate[j] = (ftr_type)malloc(sizeof(ftr_element_type) * FTR_DIM);
					}
				}
				for(r_count = 0; r_count < REPLACE_TRIAL; r_count++) {
					candidate[r_count] = random() % ds_sample->num_data;
				}
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(r_count = 0; r_count < REPLACE_TRIAL; r_count++) {
					make_pivot_center_for_QBP(piv_candidate[r_count], ds_sample->ftr_id[candidate[r_count]].ftr, med);
					piv_diff[r_count] = difference_QBP(piv_candidate[r_count], pivot->p[dim]);
				}
				piv_diff_min = FTR_DIM;
				for(r_count = 0; r_count < REPLACE_TRIAL; r_count++) {
					if(piv_diff[r_count] < piv_diff_min) {
						piv_diff_min = piv_diff[r_count];
						piv_diff_min_n = r_count;
					}
				}
				if(piv_temp == NULL) {
					piv_temp = (ftr_type)malloc(sizeof(ftr_element_type) * FTR_DIM);
				}
				// 現在のピボットを退避する．
				memcpy(piv_temp, pivot->p[dim], sizeof(ftr_element_type) * FTR_DIM);
				rad_temp = pivot->r[dim];
				memcpy(pivot->p[dim], piv_candidate[piv_diff_min_n], sizeof(ftr_element_type) * FTR_DIM);
				pivot->r[dim] = compute_rad(piv_candidate[piv_diff_min_n], dim, &sample_for_radius);  // 半径計算 (BP と同じ計算方法の場合)
			#else // REPLACE_WHOLE
				#ifdef VAR_FLIP
					num_flips = (int)(0.9 + sqrt((double)(NUM_TRIAL_AIR + NUM_TRIAL_LS - i) / (NUM_TRIAL_AIR + NUM_TRIAL_LS)) * f_max);
					if(num_flips <= 0) num_flips = 1;
				#endif
				shuffle(axi, FTR_DIM, num_flips);
				for(int j = 0; j < num_flips; j++) {
					// pivot->p[dim][axi[j]] = 255 - pivot->p[dim][axi[j]];	// 0 <=> 255 のフリップ
					pivot->p[dim][axi[j]] = (pivot->p[dim][axi[j]] == FTR_MIN ? FTR_MAX : FTR_MIN);	// FTR_MIN <=> FTR_MAX のフリップ
				}
				rad_temp = pivot->r[dim];
				pivot->r[dim] = compute_rad(pivot->p[dim], dim, &sample_for_radius);  // 半径計算 (BP と同じ計算方法の場合)
			#endif
		}
		make_y_sketch(work, dim);
		make_idx_bkt(work->y_sketch, work->ds_sample->num_data, &work->y_bucket->idx, &work->y_bucket->bkt);

//		for(int i = 0; i < resampling_size; i++) {
//			work->query_sketch[i].computed = 0;
//		}
//		make_query_sketch_for_resample(resampling_size, work);

		remake_query_sketch_for_resample(dim, resampling_size, work);
		accept = 0;
		sketch_type *temp_sketch;
		struct_bucket *temp_bucket;
//		if(obj > 0 && obj < 90) {
		temp_sketch =work->sketch; work->sketch = work->y_sketch; work->y_sketch = temp_sketch; // 精度評価のために sketch と y_sketch を交換
		temp_bucket = work->bucket; work->bucket = work->y_bucket; work->y_bucket = temp_bucket; // 精度評価のために bucket と y_bucket を交換
		eval = EVAL_PRECISION(obj, resampling_size, mode, &prec, work);
//		#if PJT_DIM < 64
//		printf("(5) Precision for K = %.1lf%% = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", (double)obj / 10, eval, prec_best, resampling_size, i);
//		#else
//		printf("(5) Precision for K = %d ppm = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", obj, eval, prec_best, resampling_size, i);
//		#endif
		if(eval >= eval_best) {
			if(final_check) {
				if(eval_best == eval) {
					not_improved++;
				} else {
					not_improved = 0;
				}
			}
			accept = 1;
			eval_best = eval;
			prec_best = prec;
		} else if(final_check) {
			not_improved++;
		}
		if(i % 100 == 0) {
			#if PJT_DIM < 64
			printf("(4) Precision for K = %.1lf%% = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", (double)obj / 10, eval_best, prec_best, resampling_size, i);
			#else
			printf("(4) Precision for K = %d ppm = %10.8lf%%, prec = %10.8lf, size = %5d, i = %6d\n", obj, eval_best, prec_best, resampling_size, i);
			#endif
			#ifdef INTERMEDIATE
				// make_query_sketch_for_resample(max_num_queries, work);
				make_query_sketch_for_resample(max_num_queries, work);
				eval = EVAL_PRECISION(obj, max_num_queries, mode, &prec, work);
				printf("max_num_queries = %d\n", max_num_queries);
				printf("i = %6d, eval = %10.8lf%%, prec = %10.8lf, INTERMEDIATE\n", i, eval, prec);
			#endif
		}
//		} else {
//			search_K_multi(multi_K);
//			if(eval_best > multi_K[obj*10]) {
//				accept = 1;
//				eval_best = multi_K[obj*10];
//			}
//		}
		if(!accept) {
//			fprintf(stderr, "Not accepted\n");
			// 元に戻す．
			#ifdef REPLACE_WHOLE
				// memcpy(pivot->p[dim], piv_temp, sizeof(ftr_element_type) * FTR_DIM);
			#else
				for(int j = 0; j < num_flips; j++) {
					pivot->p[dim][axi[j]] = (pivot->p[dim][axi[j]] == FTR_MIN ? FTR_MAX : FTR_MIN);	// FTR_MIN <=> FTR_MAX のフリップ
				}
			#endif
			pivot->r[dim] = rad_temp;
			temp_sketch = work->sketch; work->sketch = work->y_sketch; work->y_sketch = temp_sketch; // 元に戻す．
			temp_bucket = work->bucket; work->bucket = work->y_bucket; work->y_bucket = temp_bucket; // 元に戻す．
			for(int i = 0; i < resampling_size; i++) {
				work->query_sketch[i].computed = 0;
			}
			make_query_sketch_for_resample(resampling_size, work);
			remake_query_sketch_for_resample(dim, resampling_size, work);
		} else {
//			fprintf(stderr, "Accepted\n");
			num_accept++;
		}
		if(final_check && not_improved > PJT_DIM * FTR_DIM) {
			fprintf(stderr, "Convergence at i = %d\n", i - PJT_DIM * FTR_DIM);
			printf("Convergence at i = %d\n", i - PJT_DIM * FTR_DIM);
			break;
		}
	}
	
	fflush(stdout);
//	free(ans);
//	free(sketch);
//	free(y_sketch);
//	sketch = y_sketch = NULL;
	fprintf(stderr, "Num_accept = %6d\n", num_accept);
	fprintf(stderr, "Final evaluation = %10.8lf\n", eval_best);
//	printf("\nFinal evaluation = %10.8lf\n", eval_best);
}

// サンプルデータセットと質問のftrファイルと正解のcsvファイルからpivot選択のための作業域を用意する
// partition_type は，当面 QBP であるので，2を与えること．
// ftrファイルは，sample_dataset_file と query_file から，それぞれ ds_sample と ds_query に読み込む．
// サンプルデータセットから中央値medを求める．
// query_sketch を作成しておく．
// answer_file（csv）から正解情報を correct_answer に読み込む．
// query_sketch の answer に正解情報（最近傍のデータ番号と距離（sequential filteringのときはscore），answer_sketch に正解のスケッチを設定する．
// answer, sketch, y_sketch, bucket, y_bucket のメモリ確保をする．
struct_work_for_pivot_selection *new_work_for_pivot_selection(int partition_type, char *sample_dataset_file, char *query_file, char *answer_file)
{
	struct_work_for_pivot_selection *work = (struct_work_for_pivot_selection *)malloc(sizeof(struct_work_for_pivot_selection));
	work->pivot = new_pivot(partition_type);
	work->ds_sample = read_dataset_n(1, &sample_dataset_file);
	int num_data = work->ds_sample->num_data;
	fprintf(stderr, "read sample dataset OK, num_data = %d\n", num_data);
	work->ds_query = read_dataset_n(1, &query_file);
	int num_queries = work->ds_query->num_data;
	fprintf(stderr, "read query OK, num_queries = %d\n", num_queries);
	work->med = get_median(work->ds_sample);
	work->query_sketch = (struct_query_sketch *)malloc(sizeof(struct_query_sketch) * num_queries);
	for(int i = 0; i < num_queries; i++) {
		work->query_sketch[i].query = (query_type) { i, work->ds_query->ftr_id[i].ftr };
		work->query_sketch[i].computed = 0; // selected（再サンプリングで選択されたかどうか）を 0 に初期化しておく
	}
	make_query_sketch_for_resample(num_queries, work);
	work->correct_answer = read_correct_answer(answer_file, num_queries);
	for(int i = 0; i < num_queries; i++) {
		work->query_sketch[i].answer.data_num = work->correct_answer[i].data_num;
		work->query_sketch[i].answer.dist = work->correct_answer[i].dist;
		if(i < 10) printf("i = %d, data_num = %d, dist = %d\n", i, work->correct_answer[i].data_num, work->correct_answer[i].dist);
		#ifndef EXPANDED_SKETCH
		work->query_sketch[i].answer_sketch = data_to_sketch(work->query_sketch[i].query.ftr, work->pivot);
		#else
		data_to_sketch(work->query_sketch[i].query.ftr, work->pivot, work->query_sketch[i].answer_sketch);
		#endif
	}
	work->answer =  (answer_type *)malloc((num_queries * 2) * sizeof(answer_type)); // なぜ質問数の2倍？
	work->sketch = (sketch_type *)malloc(sizeof(sketch_type) * num_data);
	work->y_sketch = (sketch_type *)malloc(sizeof(sketch_type) * num_data);
//	fprintf(stderr, "sketch = %p, y_sketch = %p\n", work->sketch, work->y_sketch);
	work->bucket = (struct_bucket *)malloc(sizeof(struct_bucket));
	work->bucket->num_data = num_data;
	work->bucket->ftr_data = NULL;
	work->bucket->sk = NULL;
	work->bucket->idx = NULL;
	#ifdef NARROW_SKETCH
	work->bucket->bkt = NULL;
	#endif
	work->bucket->num_nonempty_buckets = 0;
	work->bucket->sk_num = NULL;

	work->y_bucket = (struct_bucket *)malloc(sizeof(struct_bucket));
	work->y_bucket->num_data = num_data;
	work->y_bucket->ftr_data = NULL;
	work->y_bucket->sk = NULL;
	work->y_bucket->idx = NULL;
	#ifdef NARROW_SKETCH
	work->y_bucket->bkt = NULL;
	#endif
	work->y_bucket->num_nonempty_buckets = 0;
	work->y_bucket->sk_num = NULL;

	return work;
}

// workのds_sampleのデータのスケッチをすべて求める
// for all i do sketch[i] = sketch(pivot, ftr[i]);
// narrow sketch のみ対応（wide や expanded は未対応）
void make_sketch(struct_work_for_pivot_selection *work)
{
    pivot_type *pivot = work->pivot;
    struct_dataset *ds_sample = work->ds_sample;
    sketch_type *sketch = work->sketch;
	int num_data = work->ds_sample->num_data;

	for(int i = 0; i < num_data; i++) {
		#ifndef EXPANDED_SKETCH
		sketch[i] = 0;
		#else
		for(int j = 0; j < SKETCH_SIZE; j++) {
			sketch[i][j] = 0;
		}
		#endif
	}

	#if SKETCH_PARTITION != 3
		for(int j = 0; j < PJT_DIM; j++) {
			set_dist_L2_22(pivot->p[j]);
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < num_data; i++) {
				#if defined(NARROW_SKETCH)
				write_bit(j, dist_L2_22(ds_sample->ftr_id[i].ftr) <= pivot->r[j], &sketch[i]);
				#elif defined(WIDE_SKETCH)
				write_bit_long(j, dist_L2_22(ds_sample->ftr_id[i].ftr) <= pivot->r[j], &sketch[i]);
				#else
				write_bit_expanded(j, dist_L2_22(ds_sample->ftr_id[i].ftr) <= pivot->r[j], sketch[i]);
				#endif
			}
		}
	#else
		for(int j = 0; j < PJT_DIM; j++) {
			set_part_dist_L2_22(pivot->p[j], PART_START(j), PART_DIM(j));
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i = 0; i < num_data; i++) {
				#if defined(NARROW_SKETCH)
				write_bit(j, part_dist_L2_22(ds_sample->ftr_id[i].ftr, PART_START(j), PART_DIM(j)) <= pivot->r[j], &sketch[i]);
				#elif defined(WIDE_SKETCH)
				write_bit_long(j, part_dist_L2_22(ds_sample->ftr_id[i].ftr, PART_START(j), PART_DIM(j)) <= pivot->r[j], &sketch[i]);
				#else
				write_bit_expanded(j, part_dist_L2_22(ds_sample->ftr_id[i].ftr, PART_START(j), PART_DIM(j)) <= pivot->r[j], sketch[i]);
				#endif
			}
		}
	#endif
}

// すべての i = 0, ... , num_data に対して，sketch[i] を y_sketch[i] にコピーして，dim で指定したビットだけ y_sketch[i] を書き直す．
void make_y_sketch(struct_work_for_pivot_selection *work, int dim)
{
    pivot_type *pivot = work->pivot;
    struct_dataset *ds_sample = work->ds_sample;
    sketch_type *sketch = work->sketch;
    sketch_type *y_sketch = work->y_sketch;
	int num_data = work->ds_sample->num_data;

	#if SKETCH_PARTITION != 3
		set_dist_L2_22(pivot->p[dim]);
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i = 0; i < num_data; i++) {
			#if defined(NARROW_SKETCH)
			y_sketch[i] = sketch[i];
			write_bit(dim, dist_L2_22(ds_sample->ftr_id[i].ftr) <= pivot->r[dim], &y_sketch[i]);
			#elif defined(WIDE_SKETCH)
			y_sketch[i] = sketch[i];
			write_bit_long(dim, dist_L2_22(ds_sample->ftr_id[i].ftr) <= pivot->r[dim], &y_sketch[i]);
			#else
			for(int j = 0; j < SKETCH_SIZE; j++) {
				y_sketch[i][j] = sketch[i][j];
			}
			write_bit_expanded(dim, dist_L2_22(ds_sample->ftr_id[i].ftr) <= pivot->r[dim], y_sketch[i]);
			#endif
		}
	#else
		set_part_dist_L2_22(pivot->p[dim], PART_START(dim), PART_DIM(dim));
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i = 0; i < num_data; i++) {
			#if defined(NARROW_SKETCH)
			y_sketch[i] = sketch[i];
			write_bit(dim, part_dist_L2_22(ds_sample->ftr_id[i].ftr, PART_START(dim), PART_DIM(dim)) <= pivot->r[dim], &y_sketch[i]);
			#elif defined(WIDE_SKETCH)
			y_sketch[i] = sketch[i];
			write_bit_long(dim, part_dist_L2_22(ds_sample->ftr_id[i].ftr, PART_START(dim), PART_DIM(dim)) <= pivot->r[dim], &y_sketch[i]);
			#else
			for(int j = 0; j < SKETCH_SIZE; j++) {
				y_sketch[i][j] = sketch[i][j];
			}
			write_bit_expanded(dim, part_dist_L2_22(ds_sample->ftr_id[i].ftr, PART_START(dim), PART_DIM(dim)) <= pivot->r[dim], y_sketch[i]);
			#endif
		}
	#endif
}

/*
void remake_sketch(int dim)
{
	int i;
	#pragma omp parallel for
	for(i = 0; i < db_header->data_num; i++) {
#if SKETCH_PARTITION != 3
		write_bit(dim, DISTANCE(database[i], pivot.piv0[dim], db_header->data_dim)  <= pivot.rad[dim], &sketch[i]);
#else
		write_bit(dim, pdist_L1(database[i], pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &sketch[i]);
#endif
	}
}

void make_sample_sketch(int dim, int samplesize)
{
	int i;
	#pragma omp parallel for
	for(i = 0; i < samplesize; i++)
		data_to_sketch_1bit(sample[i], dim, i);
}

void make_dbsample_sketch(int dim, int samplesize)
{
	int i;
	#pragma omp parallel for
	for(i = 0; i < samplesize; i++)
		data_to_sketch_1bit(database_temp[i], dim, i);
}

void make_query_sketch(int num_queries, struct_query_sketch qs[], pivot_type *pivot)
{
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(int i = 0; i < num_queries; i++) {
		qs[i].sketch = data_to_sketch(qs[i].query.ftr, pivot);
	}
}
*/

// work内の質問スケッチの配列query_sketchの先頭resampling_size個の設定をする．
// bd: 分割境界との最小距離, idx: bdの順位表, tbl: scoreの表計算のための表（sequential filtering のときのみ）
void make_query_sketch_for_resample(int resampling_size, struct_work_for_pivot_selection *work)
{
    pivot_type *pivot = work->pivot;
    struct_query_sketch *query_sketch = work->query_sketch;
    struct_dataset *ds_sample = work->ds_sample;
	int num_queries = work->ds_query->num_data;
//	int computed = 0;

	#pragma omp parallel for //reduction (+: computed)
	for(int i = 0; i < resampling_size; i++) {
		if(query_sketch[i].computed) continue;	// 前回の再サンプルで選択されているもののquery_sketchは計算不要
//		computed++;
		for(int p = 0; p < TABLE_SIZE; p++) {
			for(int n = 0; n < 256; n++) {
				query_sketch[i].tbl[p][n] = 0;
			}
		}
		#ifndef EXPANDED_SKETCH
		query_sketch[i].answer_sketch = 0;	// 正解（質問の最近傍）のスケッチ
		query_sketch[i].sketch = 0;			// 質問のスケッチ
		#else
		for(int j = 0; j < SKETCH_SIZE; j++) {
			query_sketch[i].answer_sketch[j] = 0;
			query_sketch[i].sketch = 0;
		}
		#endif
	}

	for(int j = 0; j < PJT_DIM; j++) {
		#if SKETCH_PARTITION != 3
		set_dist_L2_22(pivot->p[j]);
		#else
		set_part_dist_L2_22(pivot->p[j], PART_START(j), PART_DIM(j));
		#endif
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i = 0; i < resampling_size; i++) {
			if(query_sketch[i].computed) continue;	// 前回の再サンプルで計算済みのquery_sketchは計算不要

			// 質問のスケッチの第 j ビットを求める
			#if SKETCH_PARTITION != 3
			dist_type dist = dist_L2_22(query_sketch[i].query.ftr);
			#else
			dist_type dist = part_dist_L2_22(query_sketch[i].query.ftr, PART_START(j), PART_DIM(j));
			#endif
			#if defined(NARROW_SKETCH)
			write_bit(j, dist <= pivot->r[j], &(query_sketch[i].sketch));
			#elif defined(WIDE_SKETCH)
			write_bit_long(j, dist <= pivot->r[j], &(query_sketch[i].sketch));
			#else
			write_bit_expanded(j, dist <= pivot->r[j], query_sketch[i].sketch);
			#endif

			// 質問のスケッチの第 j ビットに対応する質問点から分割境界までの最短距離（＝質問とピボットの中心との距離とピボットの半径の差の絶対値）
			query_sketch[i].bd[j] = abs(dist - pivot->r[j]);
			// bd をソートする（挿入法でidxを用いた相対ソート）
			int l;
			for(l = j - 1; l >= 0 && query_sketch[i].bd[query_sketch[i].idx[l]] > query_sketch[i].bd[j]; l--) {
				query_sketch[i].idx[l + 1] = query_sketch[i].idx[l];
			}
			query_sketch[i].idx[l + 1] = j;

			// priority（score）を求めるための表関数の設定
			int tn = j / 8, bp = j % 8; // 表関数の変更（bdの追加）をすべき表番号が tn，そのビット位置が bp
			for(int n = 0; n < 256; n++) {
				if(n & (1 << bp)) {
					query_sketch[i].tbl[tn][n] += query_sketch[i].bd[j];
				}
			}

			// 正解（質問の最近傍）のスケッチの第 j ビットを求める
			#if SKETCH_PARTITION != 3
			dist_type ans_dist = dist_L2_22(ds_sample->ftr_id[query_sketch[i].answer.data_num].ftr);
			#else
			dist_type ans_dist = part_dist_L2_22(ds_sample->ftr_id[query_sketch[i].answer.data_num].ftr, PART_START(j), PART_DIM(j));
			#endif
			#if defined(NARROW_SKETCH)
			write_bit(j, ans_dist <= pivot->r[j], &(query_sketch[i].answer_sketch));
			#elif defined(WIDE_SKETCH)
			write_bit_long(j, ans_dist <= pivot->r[j], &(query_sketch[i].answer_sketch));
			#else
			write_bit_expanded(j, ans_dist <= pivot->r[j], query_sketch[i].answer_sketch);
			#endif
		}
	}
	#ifdef EVAL_BY_SEQUENTIAL_FILTERING
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(int i = 0; i < resampling_size; i++) {
		if(query_sketch[i].computed) continue;	// 前回の再サンプルで選択されているもののquery_sketchは計算不要
		query_sketch[i].computed = 1;
		query_sketch[i].answer.dist = priority(query_sketch[i].answer_sketch, &query_sketch[i]);
	}
	#endif
	
	for(int i = resampling_size; i < num_queries; i++) {
		query_sketch[i].computed = 0;
	}
//	fprintf(stderr, "resampling size = %d, computed = %d\n", resampling_size, computed); //getchar();
}

/*
void remake_query_sketch(int dim)
{
	int i;

	#pragma omp parallel for
	for(i = 0; i < num_queries; i++) {
#if SKETCH_PARTITION != 3
		write_bit(dim, DISTANCE(query[i].data, pivot.piv0[dim], db_header->data_dim)  <= pivot.rad[dim], &query[i].sketch);
#else
		write_bit(dim, pdist_L1(query[i].data, pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &query[i].sketch);
#endif
	}
}
*/

// dimで指定した次元に対応するスケッチを変更し，対応する距離下限も変更して，その順位も変更する
void remake_query_sketch_for_resample(int dim, int resampling_size, struct_work_for_pivot_selection *work)
{
    pivot_type *pivot = work->pivot;
    struct_query_sketch *query_sketch = work->query_sketch;
//	#ifdef EVAL_BY_SEQUENTIAL_FILTERING
    struct_dataset *ds_sample = work->ds_sample;
    answer_type *correct_answer = work->correct_answer;
//	#endif
	#if SKETCH_PARTITION != 3
	set_dist_L2_22(pivot->p[dim]);
	#else
	set_part_dist_L2_22(pivot->p[dim], PART_START(dim), PART_DIM(dim));
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(int i = 0; i < resampling_size; i++) {
		dist_type dist;
//		int l;

		// dim-th bit of sketch
		#if SKETCH_PARTITION != 3
		dist = dist_L2_22(query_sketch[i].query.ftr);
		#else
		dist = part_dist_L2_22(query_sketch[i].query.ftr, PART_START(dim), PART_DIM(dim));
		#endif
		#if defined(NARROW_SKETCH)
		write_bit(dim, dist <= pivot->r[dim], &(query_sketch[i].sketch));
		#elif defined(WIDE_SKETCH)
		write_bit_long(dim, dist <= pivot->r[dim], &(query_sketch[i].sketch));
		#else
		write_bit_expanded(dim, dist <= pivot->r[dim], query_sketch[i].sketch);
		#endif

		// bd
//		#ifdef EVAL_BY_SEQUENTIAL_FILTERING
			dist_type bd_old = query_sketch[i].bd[dim];
			dist_type bd_new = abs(dist - pivot->r[dim]);
			int tn = dim / 8, bp = dim % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
			for(int n = 0; n < 256; n++) {
				if(n & (1 << bp)) {
					query_sketch[i].tbl[tn][n] += (bd_new - bd_old);
				}
			}
			query_sketch[i].bd[dim] = bd_new;
			#if SKETCH_PARTITION != 3
			dist = dist_L2_22(ds_sample->ftr_id[correct_answer[query_sketch[i].query.query_num].data_num].ftr);
			#else
			dist = part_dist_L2_22(ds_sample->ftr_id[correct_answer[query_sketch[i].query.query_num].data_num].ftr, PART_START(dim), PART_DIM(dim));
			#endif
			#if defined(NARROW_SKETCH)
			write_bit(dim, dist <= pivot->r[dim], &(query_sketch[i].answer_sketch));
			#elif defined(WIDE_SKETCH)
			write_bit_long(dim, dist <= pivot->r[dim], &(query_sketch[i].answer_sketch));
			#else
			write_bit_expanded(dim, dist <= pivot->r[dim], query_sketch[i].answer_sketch);
			#endif
			query_sketch[i].answer.dist = priority(query_sketch[i].answer_sketch, &query_sketch[i]);
//		#else
//			query_sketch[i].bd[dim] = abs(dist - pivot->r[dim]);
//		#endif

		// idx (order of bd)
		for(int l = 0; l < PJT_DIM; l++) {
			if(query_sketch[i].idx[l] == dim) {
				if(l == 0 || query_sketch[i].bd[query_sketch[i].idx[l - 1]] < query_sketch[i].bd[dim]) { // 右に挿入
					for( ; l < PJT_DIM - 1 && query_sketch[i].bd[query_sketch[i].idx[l + 1]] < query_sketch[i].bd[dim]; l++) {
						query_sketch[i].idx[l] = query_sketch[i].idx[l + 1];
					}
					query_sketch[i].idx[l] = dim;
				} else { // 左に挿入
					for( ; l > 0 && query_sketch[i].bd[query_sketch[i].idx[l - 1]] > query_sketch[i].bd[dim]; l--) {
						query_sketch[i].idx[l] = query_sketch[i].idx[l - 1];
					}
					query_sketch[i].idx[l] = dim;
				}
				break;
			}
		}
	}
}

// work内のqueryからresample_size個を乱択する．
// 選んだものは質問スケッチの配列（query_sketch）の先頭にある．
void resample_query(int resample_size, struct_work_for_pivot_selection *work)
{
	int num_queries = work->ds_query->num_data;
    struct_query_sketch *query_sketch = work->query_sketch;

	for(int i = 0; i < resample_size; i++){
		int x = random() % (num_queries - i) + i;
		struct_query_sketch temp = query_sketch[i];
		query_sketch[i] = query_sketch[x];
		query_sketch[x] = temp;
	}
}

#ifdef NARROW_SKETCH
// カーソルを使った優先度付き待ち行列を使用するimamuraバージョン（delayed enq を用いている）
// 一部は，樋口のideaによる
// obj: 目的関数の指定（精度）   1 <= obj <= 69 → 候補数 k' = num_data * obj / 1000 で検索したときの精度）
//                    （候補数）90 <= obj < 100 → 精度が obj (%) になる候補数 K(%)
// ただし，現状では上の精度でのみデバッグ・調整中
// resampline_size: 評価に使用する再サンプル数（質問数）
// mode: 評価の集計方法（0: 正解数の割合，1: 正解の重み付き合計）
// *prec:
// work: 作業域
/*
double precision_resample_cursor_2_n(int obj, int resampling_size, int mode, double *prec, struct_work_for_pivot_selection *work)
{
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
//    pivot_type *pivot = work->pivot;
    struct_dataset *ds_sample = work->ds_sample;
//    struct_dataset *ds_query = work->ds_query;
    struct_query_sketch *query_sketch = work->query_sketch;
//    answer_type *correct_answer = work->correct_answer;
    answer_type *answer = work->answer;
//    sketch_type *sketch = work->sketch;
//    sketch_type *y_sketch = work->y_sketch;
    struct_bucket *bucket = work->bucket;
	int i, q;
	int num_data = ds_sample->num_data; // サンプルデータセットのデータ数
	int trunc_K;
	if(PJT_DIM < 64)
		trunc_K = (double)obj / 1000 * num_data; // 検索を優先度順の上位 trunc_K で打ち切る
	else
		trunc_K = (double)obj / 1000000 * num_data; // 検索を優先度順の上位 trunc_K で打ち切る
//	fprintf(stderr, "TRUNCATE AT %d\n", trunc_K); exit(0);
	int num_success = 0; // 正解が得られた質問数（mode = 0 のときに使用する）
	double sum = 0; // 正解の重み付き合計（mode = 1 のときに使用する）
//	int num_enumerated = 0, num_remained = 0;


	//（ここではなく，呼び出し元のAIRでsketchの変更があったときのみ作り直すようにした方がよいかも
	// make_idx_bkt(sketch, num_data, &bucket->idx, &bucket->bkt); // 評価の度にバケット表を作り直す）
	#ifdef _OPENMP
	int nt = omp_get_max_threads();
	#pragma omp parallel for
	#else
	int nt = 1;
	#endif
	for(i = 0; i < resampling_size; i++)
		answer[i].dist = UINT_MAX;

	static struct_que_c2_n *que_pool = NULL;
	int used_q[nt];
	if(que_pool == NULL) {
		que_pool = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n) * nt);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
	}
	for(q = 0; q < nt; q++) {
		used_q[q] = 0;
	}

	int *bkt = bucket->bkt, *idx = bucket->idx;
#ifdef TRAINING_PARALLEL
	#ifdef _OPENMP
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
	#endif
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		QUE_c2 qu, qu2;
		struct_que_c2_n *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均
//		int enumerated = 0;

//	if(lb == 1 || lb >= 5) { // L1
		k = 0;
//#if SKETCH_PARTITION == 3
//		pruned = 0;
//#endif

		s = query_sketch[q].sketch;
		mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
		for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
			if(idx[j] == query_sketch[q].answer.data_num) {
//				printf("(1) found: q = %d, data_num = %d\n", q, idx[j]);
				num_success++;
				sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
				answer[q].data_num = q + 1;
				answer[q].dist = k + 1;
				goto END_SEARCH; // 検索終了
			} else if(idx[j] > query_sketch[q].answer.data_num) {
				k += bkt[s + 1] - j; 
				break;
			}
		}

		s = s ^ (1 << query_sketch[q].idx[0]);
		if(query_sketch[q].bd[query_sketch[q].idx[0]] > score) {
			score = query_sketch[q].bd[query_sketch[q].idx[0]];
			if(k >= trunc_K) goto END_SEARCH;
		}
		mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
		for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
			if(idx[j] == query_sketch[q].answer.data_num) {
//				printf("(2) found: q = %d, data_num = %d\n", q, idx[j]);
				num_success++;
				sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
				answer[q].data_num = q + 1;
				answer[q].dist = k + 1;
				goto END_SEARCH; //k = db_header->data_num; // 検索終了
			} else if(idx[j] > query_sketch[q].answer.data_num) {
				k += bkt[s + 1] - j;
				break;
			}
		}

		int tn;
		#ifdef _OPENMP
		#pragma omp critical 
		#endif
		{
			for(tn = 0; tn < nt; tn++) {
				if(used_q[tn] == 0) {
					que = &que_pool[tn];
					used_q[tn] = 1;
					break;
				}
				if(tn == nt) {
					fprintf(stderr, "No free que in que_pool\n");
					exit(0);
				}
			}
		}

		make_empty_que_c2_n(que);

		// enq pattern of 0...10
		qu.cursor = new_que_e2_n(que);
		qu.key = query_sketch[q].bd[query_sketch[q].idx[1]];
		que->details[qu.cursor].sk = query_sketch[q].sketch ^ (1 << query_sketch[q].idx[1]);
		que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
		enq_c2_n(&qu, que);

		while(deq_c2_n(&qu, que)) {
//	START_SEARCH:
//				enumerated++;
//				printf("cursor = %d, pt = %d, key = %d\n", qu.cursor, que->details[qu.cursor].pt, qu.key);
			s = que->details[qu.cursor].sk;
			if(qu.key > score) {
				if(k >= trunc_K) goto FREE_Q;
				score = qu.key;
			}
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query_sketch[q].answer.data_num) {
//					printf("(3) found: q = %d, data_num = %d\n", q, idx[j]);
					num_success++;
					sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
					answer[q].data_num = q + 1;
					answer[q].dist = k + 1;
					goto FREE_Q; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query_sketch[q].answer.data_num) {
					k += bkt[s + 1] - j;
					break;
				}
			}
			switch(que->details[qu.cursor].pt & 15) {
			case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
			case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
				{
					int m = lsb_pos(que->details[qu.cursor].pt);
					if(m > 0 && m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
						// Y010^m -> Y10^{m+1}
						qu2.cursor = new_que_e2_n(que);
						qu2.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[m + 1]] - query_sketch[q].bd[query_sketch[q].idx[m]];
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[m + 1])) ^ (1 << query_sketch[q].idx[m]);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
						// Y010^m -> Y010^{m-1}1
						qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
						enq_c2_n(&qu2, que);
					} else {
						qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
					}
				}
				break;
			case 4:  // X0100 -> enq(X0101) and enq(X1000)
				// X1000
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[3]] - query_sketch[q].bd[query_sketch[q].idx[2]];
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[3])) ^ (1 << query_sketch[q].idx[2]);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
				// X0101
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 1:  // X0001 -> enq(X0010)
			case 5:  // X0101 -> enq(X0110)
			case 9:  // X1001 -> enq(X1010)
			case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[1]] - query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[1])) ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 2:  // X0010 -> enq(X0011) and enq(X0100)
			case 10: // X1010 -> enq(X1011) and enq(X1100)
				// X0100 and X1100
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[2]] - query_sketch[q].bd[query_sketch[q].idx[1]];
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[2])) ^ (1 << query_sketch[q].idx[1]);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
				// X0011 and X1011
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 6:  // X0110 -> enq(X0111)
			case 12: // X1100 -> enq(X1101)
			case 14: // X1110 -> enq(10111)
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 3:  // X0011
			case 7:  // X0111
			case 11: // X1011
			case 15: // X1111 -> nothing to do
				break;
			}
		}
FREE_Q:
//fprintf(stderr, "FREE_Q, q = %d\n", q);
		#ifdef _OPENMP
		#pragma omp critical 
		#endif
		{
			used_q[tn] = 0;
//				num_enumerated += enumerated;
//				num_remained += que->qsize;
		}
END_SEARCH: ;
//fprintf(stderr, "END_SEARCH, q = %d\n", q);
//fprintf(stderr, "*prec = %.2lf, eval = %.2lf\n", (double)100 * num_success / resampling_size, sum / resampling_size);
	}
//fprintf(stderr, "prec = %.2lf, eval = %.2lf, mode = %d\n", (double)100 * num_success / resampling_size, sum / resampling_size, mode);
//	printf("max detail size = %d\n", max_detail_size); // exit(0);
//	printf("enumerated = %d, remained = %d, shurinked = %6.2lf%%\n", num_enumerated, num_remained, (double)(num_enumerated - num_remained) / num_enumerated); exit(0);
//	printf("num_success = %d, resampling_size = %d\n", num_success, resampling_size);

	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}
*/
double precision_resample_cursor_2_n(int obj, int resampling_size, int mode, double *prec, struct_work_for_pivot_selection *work)
{
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
    struct_dataset *ds_sample = work->ds_sample;
    struct_query_sketch *query_sketch = work->query_sketch;
    answer_type *answer = work->answer;
    struct_bucket *bucket = work->bucket;
	int i, q;
	int num_data = ds_sample->num_data; // サンプルデータセットのデータ数
	int trunc_K;
	if(PJT_DIM < 64)
		trunc_K = (double)obj / 1000 * num_data; // 検索を優先度順の上位 trunc_K で打ち切る
	else
		trunc_K = (double)obj / 1000000 * num_data; // 検索を優先度順の上位 trunc_K で打ち切る
	int num_success = 0; // 正解が得られた質問数（mode = 0 のときに使用する）
	double sum = 0; // 正解の重み付き合計（mode = 1 のときに使用する）

	#ifdef _OPENMP
	int nt = omp_get_max_threads();
	#pragma omp parallel for
	#else
	int nt = 1;
	#endif
	for(i = 0; i < resampling_size; i++)
		answer[i].dist = UINT_MAX;

	static struct_que_c2_n *que_pool = NULL;
	int used_q[nt];
	if(que_pool == NULL) {
		que_pool = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n) * nt);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
	}
	for(q = 0; q < nt; q++) {
		used_q[q] = 0;
	}

	int *bkt = bucket->bkt, *idx = bucket->idx;
#ifdef TRAINING_PARALLEL
	#ifdef _OPENMP
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
	#endif
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		QUE_c2 qu, qu2;
		struct_que_c2_n *que = NULL;
		int k;
//		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均
		k = 0;

		s = query_sketch[q].sketch;
		if(s == query_sketch[q].answer_sketch) {
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2;
			num_success++;
			sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
			answer[q].data_num = q + 1;
			answer[q].dist = k + 1;
//			printf("(1) found: q = %d, k = %d\n", q, k);
			goto END_SEARCH; // 検索終了
		}
		k += bkt[s + 1] - bkt[s];
		if(k >= trunc_K) {
//			printf("(1) not found: q = %d, k = %d\n", q, k);
			goto END_SEARCH;
		}

		s = s ^ (1 << query_sketch[q].idx[0]);
//		if(query_sketch[q].bd[query_sketch[q].idx[0]] > query_sketch[q].answer.dist) goto END_SEARCH; // s の priority が正解の priority を超えた => fail
		if(s == query_sketch[q].answer_sketch) {
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2;
			num_success++;
			sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
			answer[q].data_num = q + 1;
			answer[q].dist = k + 1;
//			printf("(2) found: q = %d, k = %d\n", q, k);
			goto END_SEARCH; // 検索終了
		}
		k += bkt[s + 1] - bkt[s];
		if(k >= trunc_K) {
//			printf("(2) not found: q = %d, k = %d\n", q, k);
			goto END_SEARCH;
		}

		int tn;
		#ifdef _OPENMP
		#pragma omp critical 
		#endif
		{
			for(tn = 0; tn < nt; tn++) {
				if(used_q[tn] == 0) {
					que = &que_pool[tn];
					used_q[tn] = 1;
					break;
				}
				if(tn == nt) {
					fprintf(stderr, "No free que in que_pool\n");
					exit(0);
				}
			}
		}

		make_empty_que_c2_n(que);

		// enq pattern of 0...10
		qu.cursor = new_que_e2_n(que);
		qu.key = query_sketch[q].bd[query_sketch[q].idx[1]];
		que->details[qu.cursor].sk = query_sketch[q].sketch ^ (1 << query_sketch[q].idx[1]);
		que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
		enq_c2_n(&qu, que);

		while(deq_c2_n(&qu, que)) {
			s = que->details[qu.cursor].sk;
//			if(qu.key > query_sketch[q].answer.dist) goto FREE_Q; // s の priority が正解の priority を超えた => fail
			if(s == query_sketch[q].answer_sketch) {
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2;
				num_success++;
				sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
				answer[q].data_num = q + 1;
				answer[q].dist = k + 1;
//				printf("(3) found: q = %d, k = %d\n", q, k);
				goto FREE_Q; // 検索終了
			}
			k += bkt[s + 1] - bkt[s];
			if(k >= trunc_K) {
//				printf("(3) not found: q = %d, k = %d", q, k);
				goto FREE_Q;
			}

			switch(que->details[qu.cursor].pt & 15) {
			case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
			case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
				{
					int m = lsb_pos(que->details[qu.cursor].pt);
					if(m > 0 && m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
						// Y010^m -> Y10^{m+1}
						qu2.cursor = new_que_e2_n(que);
						qu2.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[m + 1]] - query_sketch[q].bd[query_sketch[q].idx[m]];
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[m + 1])) ^ (1 << query_sketch[q].idx[m]);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
						// Y010^m -> Y010^{m-1}1
						qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
						enq_c2_n(&qu2, que);
					} else {
						qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
					}
				}
				break;
			case 4:  // X0100 -> enq(X0101) and enq(X1000)
				// X1000
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[3]] - query_sketch[q].bd[query_sketch[q].idx[2]];
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[3])) ^ (1 << query_sketch[q].idx[2]);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
				// X0101
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 1:  // X0001 -> enq(X0010)
			case 5:  // X0101 -> enq(X0110)
			case 9:  // X1001 -> enq(X1010)
			case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[1]] - query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[1])) ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 2:  // X0010 -> enq(X0011) and enq(X0100)
			case 10: // X1010 -> enq(X1011) and enq(X1100)
				// X0100 and X1100
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[2]] - query_sketch[q].bd[query_sketch[q].idx[1]];
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[2])) ^ (1 << query_sketch[q].idx[1]);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
				// X0011 and X1011
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 6:  // X0110 -> enq(X0111)
			case 12: // X1100 -> enq(X1101)
			case 14: // X1110 -> enq(10111)
				qu.key = qu.key + query_sketch[q].bd[query_sketch[q].idx[0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query_sketch[q].idx[0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 3:  // X0011
			case 7:  // X0111
			case 11: // X1011
			case 15: // X1111 -> nothing to do
				break;
			}
		}
FREE_Q:
		#ifdef _OPENMP
		#pragma omp critical 
		#endif
		{
			used_q[tn] = 0;
		}
END_SEARCH: ;
	}

	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}
#endif

#ifdef EVAL_BY_SEQUENTIAL_FILTERING
// AIRやLSによるピボット探索時に用いる Sequential Filtering に基づく，精度（recall）の評価		
double precision_resample_by_sequential_filtering(int obj, int resampling_size, int mode, double *prec, struct_work_for_pivot_selection *work)
{
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合（未実装）
	//                         1 → 重み付きスコア合計
    // pivot_type *pivot = work->pivot;
    struct_dataset *ds_sample = work->ds_sample;
    // struct_dataset *ds_query = work->ds_query;
    struct_query_sketch *query_sketch = work->query_sketch;
    // answer_type *correct_answer = work->correct_answer;
    answer_type *answer = work->answer;
    sketch_type *sketch = work->sketch;
    // sketch_type *y_sketch = work->y_sketch;
    // struct_bucket *bucket = work->bucket;
	int q;
	int num_data = ds_sample->num_data; // サンプルデータセットのデータ数
	int trunc_K;
	if(PJT_DIM < 64)
		trunc_K = (double)obj / 1000 * num_data; // 検索を優先度順の上位 trunc_K で打ち切る
	else
		trunc_K = (double)obj / 1000000 * num_data; // 検索を優先度順の上位 trunc_K で打ち切る
	int num_success = 0; // 正解が得られた質問数（mode = 0 のときに使用する）
	double sum = 0; // 正解の重み付き合計（mode = 1 のときに使用する）

	#pragma omp parallel for
	for(q = 0; q < resampling_size; q++)
		answer[q].dist = UINT_MAX;

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: sum, num_success) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		//	各射影次元 dim について，（dim = 0, ... , PJT_DIM - 1）分割境界と質問の最小距離を求める．
		//	bd[dim] = | dist(q, p[dim]) - r[dim] | （score_2のときは，それを2乗しておく）
		//	8ビットごとに分割した優先順位（score_p）を求めるための表関数を作成する．
		//	以上は，すでに作成済みであるとする（make_query_sketch_for_resample および remake_query_sketch_for_resample）

		// local variables for parallel
		dist_type score;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均
		int n_low = 0, n_eq = 0; // 優先順位が正解のものより高い（スコアが低い）ものの個数，正解のものと同じ（スコアが等しい）ものの個数
		int i = 0;
		for( ; i < num_data / 8; i++) {
			score = priority(sketch[i], &query_sketch[q]);
			if(score < query_sketch[q].answer.dist) {
				n_low++;
			} else if(score == query_sketch[q].answer.dist) {
				n_eq++;
			}
		}
		if(n_low > trunc_K) continue;

		for( ; i < num_data / 4; i++) {
			score = priority(sketch[i], &query_sketch[q]);
			if(score < query_sketch[q].answer.dist) {
				n_low++;
			} else if(score == query_sketch[q].answer.dist) {
				n_eq++;
			}
		}
		if(n_low > trunc_K) continue;

		for( ; i < num_data / 2; i++) {
			score = priority(sketch[i], &query_sketch[q]);
			if(score < query_sketch[q].answer.dist) {
				n_low++;
			} else if(score == query_sketch[q].answer.dist) {
				n_eq++;
			}
		}
		if(n_low > trunc_K) continue;

		for( ; i < num_data; i++) {
			score = priority(sketch[i], &query_sketch[q]);
			if(score < query_sketch[q].answer.dist) {
				n_low++;
			} else if(score == query_sketch[q].answer.dist) {
				n_eq++;
			}
		}
		if(n_low > trunc_K) continue;

		mid_k = n_low + n_eq / 2;
		if(mid_k < trunc_K) {
//			printf("q = %d, n_low = %d, n_eq = %d, answer.data_num = %d\n", q, n_low, n_eq, query_sketch[q].answer.data_num);
			num_success++;
		}
		sum += log(SCORE_FACTOR * ds_sample->num_data / (mid_k + 1));
		answer[q].data_num = q + 1;
		answer[q].dist = mid_k + 1;
	}

//	printf("num_success = %d, resampling_size = %d\n", num_success, resampling_size);

	*prec = (double)100 * num_success / resampling_size;

	return (sum / resampling_size);
}
#endif

//  n 個の整数（0, 1, ... , n - 1) から m 個を乱択するためのシャッフル
// 使い方：
// 1. 準備
//    n 個の整数（0, 1, ... , n - 1) が入った配列を用意する．（arrayとする）
// 2. shuffle(array, m)で array の要素をシャッフルする．
// 3. array[0], array[1], ... , array[m - 1] が求める整数

void shuffle(int a[], int n, int m)
{
	int i, r, t;
	for(i = 0; i < m; i++) {
		r = random() % (n - i);
//		r = random(n - i);
		t = a[i];
		a[i] = a[r + i];
		a[r + i] = t;
	}
}

#ifdef NARROW_SKETCH
// スケッチからバケット表（idxとbkt）を作成する
// sk = スケッチ，num_data = データ数
// idx, bkt = バケット表の配列へのポインタ
// スケッチ sk[0], sk[1], ... , sk[num_data - 1] をスケッチを整数と見て昇順に並べたとき（実際にソートするわけでなはい）
// bkt[s] = スケッチ s の最左位置
// idx[i] = スケッチの昇順で i 番目のスケッチの元の位置（スケッチを相対ソートしたときのインデックスに相当する）
// つまり，idx[0], ... , idx[num_data - 1] は，0, ... , num_data - 1 の並べ替えで，
// sk[idx[0]] <= sk[idx[1]] <= ... <= sk[idx[num_data - 1]] になっている（相対ソート）
// スケッチ s は，スケッチ順の並びで，bkt[s], bkt[s] + 1, ... , bkt[s + 1] - 1 の位置にある．
// スケッチ s は，元の並びでは，idx[bkt[s]], idx[bkt[s] + 1], ... , idx[bkt[s + 1] - 1] に現れている．
// スケッチ s の出現数は，bkt[s + 1] - bkt[s] である．bkt[s] = bkt[s + 1] のときは，s の出現数は 0 である．
void make_idx_bkt(sketch_type sk[], int num_data, int **idx, int **bkt)
{
	unsigned int bkt_size = (1U << PJT_DIM);

//	fprintf(stderr, "num_data = %d, &sk[0] = %p, *idx = %p, *bkt = %p\n", num_data, &sk[0], *idx, *bkt);

	if(*idx == NULL) *idx = (int *)malloc(sizeof(int) * num_data);
	if(*bkt == NULL) *bkt = (int *)malloc(sizeof(int) * (bkt_size + 2));

	int *bucket = *bkt + 1;
	for(int i = 0; i < bkt_size + 1; i++) bucket[i] = 0;
//	fprintf(stderr, "clear bkt OK\n");
	for(int i = 0; i < num_data; i++) bucket[sk[i] + 1]++;
//	fprintf(stderr, "count sketch OK\n");
	for(int i = 0; i < bkt_size; i++) bucket[i + 1] += bucket[i];
//	fprintf(stderr, "compute bkt OK\n");
//	for(int i = 0; i < num_data; i++) {
//		if(i < 5) {
//			fprintf(stderr, "i = %d\n", i);
//			fprintf(stderr, "sk[%d] = %d\n", i, sk[i]);
//			fprintf(stderr, "bucket[%d] = %d\n", sk[i], bucket[sk[i]]);
//		}
//		(*idx)[bucket[sk[i]]++] = i;
//	}
	for(int i = 0; i < num_data; i++) (*idx)[bucket[sk[i]]++] = i;
//	fprintf(stderr, "compute idx OK\n");
	(*bkt)[0] = 0;
//	fprintf(stderr, "*idx = %p, *bkt = %p\n", *idx, *bkt);
}
#endif

void make_pivot_center_for_QBP(ftr_type quantized, ftr_type center, ftr_type med)
{
//	quantized = center を量子化した点
//	center = データベースから乱択した点
//	med = 中央値
//	QBP のみ（BP や PQBP は未対応？)
//	int dim_from = 0, dim_to = 0;
//#if SKETCH_PARTITION == 3
//	dim_from = dim * 4;
//	dim_to = (dim + 1) * 4 - 1;
//#endif

//	if(pivot.partition == 1)   { // BP
//		for(int i = 0; i < db_header->data_dim; i++) {
//			pivot.piv0[dim][i] = database[c][i];
//		}
//		memcpy(center, database[c], db_header->element_size * db_header->data_dim);
//	} else if(pivot.partition == 2)   { // QBP
		for(int i = 0; i < FTR_DIM; i++) {
			quantized[i] = (center[i] < med[i]) ? FTR_MIN : FTR_MAX;
		}
//	} else if(pivot.partition == 3)   { // PQBP
//		for(int i = dim_from; i <= dim_to; i++) {
//			center[i] = (database[c][i] < med[i]) ? 0 : 255;
//		}
//	}
}
/*
static int comp(const void *a, const void *b)
{
	if(*((unsigned int *) a) < *((unsigned int *) b))
		return -1;
	else if(*((unsigned int *) a) == *((unsigned int *) b))
		return 0;
	else
		return 1;
}

dist_type compute_rad(ftr_type center, int dim, int samplesize, int *count_tbl)
{
	int i;
	dist_type result;

	for(i = 0; i < samplesize; i++){	
		count_tbl[i] = DISTANCE(sample[i], piv, dim);
	}
	qsort(count_tbl, samplesize, sizeof(dist_type), comp);
	int x;
	if(samplesize % 2 == 1) {
		x = (samplesize + 1) / 2;
		result = count_tbl[x];
	}
	else {
		x = samplesize / 2;
		result = (count_tbl[x] + count_tbl[x + 1]) / 2;
	}
	return result;
}
*/

// ピボットの半径を求める．
// center 		= ピボットの中心点
// dim 			= ピボットの次元番号（スケッチの何番目のビットかを表すもので，スケッチ幅 w ではない） 
//                PQBP のときに部分距離を求めるために使用する．
// sample_size 	= 半径を求めるためのサンプル点の個数
// sample 		= サンプル点
// 【方法】center と sample[0], ... , sample[sample_size - 1] の距離の中央値を求める．
dist_type compute_rad(ftr_type center, int dim, ftr_sample *sample)
{
	dist_type r;
	static dist_type *sample_dist = NULL;
	if(sample_dist == NULL) {
		sample_dist = (dist_type *)malloc(sizeof(dist_type) * sample->num_samples);
	}
	static int *idx = NULL;
	if(idx == NULL) {
		idx = (int *)malloc(sizeof(int) * sample->num_samples);
	}

	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(int i = 0; i < sample->num_samples; i++) {
		#ifndef PARTITION_TYPE_PQBP
		sample_dist[i] = DISTANCE(center, sample->ftr[i], FTR_DIM);
		#else
		sample_dist[i] = PART_DISTANCE(center, sample->ftr[i], PART_START(dim), PART_DIM(dim));
		#endif
		idx[i] = i;
	}
	#ifdef _OPENMP
	int nt = omp_get_max_threads();
	#else
	int nt = 1;
	#endif
	if(sample->num_samples >= 1000) {
		r = get_threshold_k(idx, sample_dist, sample->num_samples, sample->num_samples / 2, nt);  // 半径計算(get_thresholdを用いて中央値を求めている)
//		r = get_threshold_k(idx, sample_dist, sample->num_samples, 3 * sample->num_samples / 9 + random() % (2 * sample->num_samples / 9), nt);  // 半径計算(中央値より少しずらしたものを半径にする)
	} else {
		insertion_sort(sample_dist, sample->num_samples); // 半径計算に使用するサンプル数が少ないときは insertion_sort を用いて中央値を計算する
		r = sample_dist[sample->num_samples / 2];
	}

	return r;
}

int difference_QBP(ftr_type c1, ftr_type c2) {
	int i, d = 0;
	for(i = 0; i < FTR_DIM; i++) {
		d += (c1[i] != c2[i]);
	}
	return d;
}
#endif // USE_AIR