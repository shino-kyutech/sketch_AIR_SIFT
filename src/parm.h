#ifndef COMPILE_TIME // 実際のコンパイル時にはCOMPILE_TIMEをdefineして，以下の定義を無効にすること．
                     // VSCODEでプログラムを作成している時に，未定義定数でエラーが出るのを防ぐために，
                     // 実際のコンパイル時に事前に定義して用いる定数を定義しておく． 

#define DEEP1B
#define NUM_K 10    // k-NN で求める近傍数
#define QUERY "query_xx.ftr"    // 質問データのftrファイル
#define SAMPLE_DATASET "base10M_00.ftr"     // ピボット選択時のサンプルデータセット
#define PIVOT_FILE "pivot.csv"
#define BUCKET_FILE "pivot_range.bkt"
#define QUERY_FILE "query_xx.ftr"
#define ANSWER_FILE "query_xx.csv"
#define QUERY_2ND_FILE "query_xx.ftr"
#define ANSWER_2ND_FILE "query_xx.csv"
#define QUERY_3RD_FILE "query_xx.ftr"
#define ANSWER_3RD_FILE "query_xx.csv"
#define SFTR_FILE "base_yy_xx.ftr"
#define SEED 1 // random seed

#define PARTITION_TYPE_QBP
//#define PARTITION_TYPE_PQBP

#define _OPENMP
#define NUM_THREADS 8
#define NARROW_SKETCH
//#define WIDE_SKETCH
//#define EXPANDED_SKETCH

#define PJT_DIM 16
#define SAMPLE_SIZE 10000
//#define FTR_ON_MAIN_MEMORY
#define FTR_ON_SECONDARY_MEMORY

#define FTR_ARRANGEMENT_ASIS

#define RESULT_FILE "result/${pivot}_${range}.csv"
#define NUM_CANDIDATES  100 // データセットに対する割合（パーミリアド（万分率））: 100 -> 1%
#define NUM_CANDIDATES1 200 // データセットに対する割合（パーミリアド（万分率））: 200 -> 2%
#define NUM_CANDIDATES2 300
#define NUM_CANDIDATES3 400
#define NUM_CANDIDATES4 500
#define NUM_CANDIDATES5 600

//#define FILTERING_BY_SKETCH_ENUMERATION
#define FILTERING_BY_SKETCH_ENUMERATION_C2N
//#define SEQUENTIAL_FILTERING
//#define SEQUENTIAL_FILTERING_USING_BUCKET
//#define SEQUENTIAL_FILTERING_USING_HAMMING
//#define DOUBLE_FILTERING

#define SELECT_K_BY_QUICK_SELECT_K_PARA_WORK

#define NUM_Q 100
#define NUM_D 10000
#define REPEAT 100
#define ACCESS_MODE "SEQUENTIAL"

#define EXPANDED_PIVOT_FILE     "pivot_PQBP_t10_w192_sd00_ss10000_np6_sr1000_seed1.csv"
#define EXPANDED_BUCKET_FILE    "pivot_PQBP_t10_w192_sd00_ss10000_np6_sr1000_seed1_00_99.bkt"
#define NARROW_PIVOT_FILE       "pivot_QBP_t1000_w28_sd00_ss10000.csv"
#define NARROW_BUCKET_FILE      "pivot_QBP_t1000_w28_sd00_ss10000_00_99.bkt"

#define NUM_CANDIDATES_1ST      3000
#define NUM_CANDIDATES_1ST_END  8000
#define NUM_CANDIDATES_1ST_STEP 1000
#define NUM_CANDIDATES_2ND      1
#define NUM_CANDIDATES_2ND_END  5
#define NUM_CANDIDATES_2ND_STEP 1

#define SCORE_P_1ST 1.0
#define SCORE_P_2ND 1.0

#define PRIORITY 0

//#define SCORE_1
//#define PRIORITY 1

//#define SCORE_2
//#define PRIORITY 2

//#define SCORE_INF
//#define PRIORITY 3

#define NUM_TRIAL_QBP 100 // Number of Trials in incremental QBP
// Macro for AIR
#define USE_AIR
#ifdef USE_AIR
#define NUM_TRIAL_AIR 500 // Number of Trials in AIR
#define TRUNCATE_AT 10    // Number of candidates (% to num_data) to use evaluate pivots in optimization by AIR 
#define NUM_TRIAL_LS  0	  // Number of Trials in Local Search
#define MAX_SAMPLESIZE 10000
#define RETRY 200
#define REUSE 10	//取り換え頻度
#define T_POWER 2
#define N0 500				//初期サンプルサイズ

#define EVAL_BY_SEQUENTIAL_FILTERING // AIRでピボットの評価に sequential filtering を用いるときに #define

#define EVAL_PRECISION precision_resample_cursor_2_n
// #define EVAL_PRECISION precision_resample_by_sequential_filtering

#define K1 0.6645		    //95%到達…500-0.2802 100-0.5488
//#define K1 0.2802
#define K2 1.0
// #define K 131028

// #define SKIP_ENQ		// score による打ち切りを行うときに #define 
//#ifdef SKIP_ENQ
//double expand = 1.0 ; // 射影距離を引き延ばしたり，縮めたりする係数
//#define SKIP_CONDITION(pdist, dist) ((pdist) * expand > (dist))
// スケッチの列挙で用いる優先度付き待ち行列にスケッチをENQ（挿入）するときに，
// 射影距離 pdist × expand がその時点での最近傍暫定解との距離 dist を超えていたら，ENQ をしない．
// DISA の DeCAF記述子 においては，実距離に比べて非常に大きいので，expand は 0.1 のように小さくすること．
//#endif

#define SCORE_FACTOR 0.4

//#define REPLACE_WHOLE
#define REPLACE_TRIAL 10

#define VAR_FLIP
#define NUM_FLIPS 1

#ifndef EVAL_MODE		// 評価にlogを使った重み付きスコア総和を用いるときは 1 にする．
#define EVAL_MODE 0
#endif
#endif
#endif
// これより下には追加しないこと