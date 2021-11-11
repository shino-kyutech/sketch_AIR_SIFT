#include "sketch.h"

#ifdef USE_AIR
//    pivot_type *pivot;                      // ピボット
//    struct_dataset *ds_sample;              // ピボット評価（検索精度（recall））のための（小さな）サンプルデータセット
//    struct_dataset *ds_query;               // 質問データセット
//    ftr_type med;                           // サンプルデータセットの中央値
//    struct_query_sketch *query_sketch;      // 検索のための構造体（質問，bd, idx, tbl, answer）
//    answer_type *correct_answer;            // 正解情報（data_num，dist）
//    answer_type *answer;                    // 検索時の解（data_num，dist または score）
//    sketch_type *sketch;                    // サンプルデータのスケッチ
//    sketch_type *y_sketch;                  // サンプルデータのスケッチ（一時退避用）
//    struct_bucket *bucket;                  // サンプルデータセットのバケット表
//    struct_bucket *y_bucket;                // サンプルデータセットのバケット表（一時退避用）
typedef struct {
    pivot_type *pivot;                      // ピボット
    struct_dataset *ds_sample;              // ピボット評価（検索精度（recall））のための（小さな）サンプルデータセット
    struct_dataset *ds_query;               // 質問データセット
    ftr_type med;                           // サンプルデータセットの中央値
    struct_query_sketch *query_sketch;      // 検索のための構造体（質問，bd, idx, tbl, answer）
    answer_type *correct_answer;            // 正解情報（data_num，dist）
    answer_type *answer;                    // 検索時の解（data_num，dist または score）
    sketch_type *sketch;                    // サンプルデータのスケッチ
    sketch_type *y_sketch;                  // サンプルデータのスケッチ（一時退避用）
    struct_bucket *bucket;                  // サンプルデータセットのバケット表
    struct_bucket *y_bucket;                // サンプルデータセットのバケット表（一時退避用）
} struct_work_for_pivot_selection;
#endif
    // pivot_type *pivot = work->pivot;
    // struct_dataset *ds_sample = work->ds_sample;
    // struct_dataset *ds_query = work->ds_query;
    // ftr_type med = work->med;
    // struct_query_sketch *query_sketch = work->query_sketch;
    // answer_type *correct_answer = work->correct_answer;
    // answer_type *answer = work->answer;
    // sketch_type *sketch = work->sketch;
    // sketch_type *y_sketch = work->y_sketch;
    // struct_bucket *bucket = work->bucket;

// int num_samples;
// ftr_type *ftr;
typedef struct {
    int num_samples;
    ftr_type *ftr;
} ftr_sample;

void get_sample(struct_dataset *ds, ftr_sample *sample);
ftr_type get_median(struct_dataset *ds);
void select_pivot_QBP(pivot_type *pivot, ftr_type median, ftr_sample *sample, struct_dataset *ds, int nt);
void select_pivot_random_QBP(pivot_type *pivot, ftr_type median, struct_dataset *ds);
double collision(int num_sk, sketch_type sk[]);
// wide スケッチ用は，まずは，QBP のみで 
// optimize_pivot_by_precision // LS and AIR
#ifdef USE_AIR
void optimize_pivot_by_precision_AIR_flip(struct_work_for_pivot_selection *work);
struct_work_for_pivot_selection *new_work_for_pivot_selection(int patrition_type, char *sample_dataset_file, char *query_file, char *answer_file);
void make_sketch(struct_work_for_pivot_selection *work);
void remake_sketch(int dim);
void make_y_sketch(struct_work_for_pivot_selection *work, int dim);
void make_sample_sketch(int dim, int samplesize);
void make_dbsample_sketch(int dim, int samplesize);
// void make_query_sketch(int num_queries, struct_query_sketch qs[]);

void make_query_sketch_for_resample(int resampling_size, struct_work_for_pivot_selection *work);
void remake_query_sketch(int dim);
void remake_query_sketch_for_resample(int dim, int resampling_size, struct_work_for_pivot_selection *work);
void resample_query(int samplesize, struct_work_for_pivot_selection *work);
double precision_resample_cursor_2_n(int num_K, int resample_size, int mode, double *prec, struct_work_for_pivot_selection *work);
#ifdef EVAL_BY_SEQUENTIAL_FILTERING
double precision_resample_by_sequential_filtering(int obj, int resampling_size, int mode, double *prec, struct_work_for_pivot_selection *work);
#endif
void shuffle(int a[], int n, int m);
void make_idx_bkt(sketch_type sk[], int num_data, int **idx, int **bkt);
void make_pivot_center_for_QBP(ftr_type quantized, ftr_type center, ftr_type med);
dist_type compute_rad(ftr_type center, int dim, ftr_sample *sample);
int difference_QBP(ftr_type c1, ftr_type c2);
#endif
