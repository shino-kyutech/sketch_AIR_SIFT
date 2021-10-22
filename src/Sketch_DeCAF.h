#pragma once

#include "Config.h"
#include "ftr.h"

#ifndef SCORE_FACTOR
#define SCORE_FACTOR 0.4
#endif

#ifndef ENQ_P
#define ENQ_P enq_p
#define DEQ_P deq_p
#endif
void print_bin(unsigned int s);
int msb_pos(unsigned int x);
int lsb_pos(unsigned int x);
int msb_pos_long(unsigned long x);
int lsb_pos_long(unsigned long x);
void median(void);
int difference_QBP(ftr_type c1, ftr_type c2);
void make_pivot_for_QBP(int candidate, ftr_type center, int dim);
void select_pivot_BPorQBP(void);
void random_QBP(void);
void optimize_pivot_by_precision_AIR_flip(int obj, int mode, int f_num);
void make_sketch(void);
void remake_sketch(int dim);
void make_y_sketch(int dim);
void make_sample_sketch(int dim, int samplesize);
void make_dbsample_sketch(int dim, int samplesize);
void make_query_sketch(void);
void make_query_sketch_for_resample(int resampling_size);
void remake_query_sketch(int dim);
void remake_query_sketch_for_resample(int dim, int resampling_size);
unsigned int data_to_sketch(ftr_type o, int dim);
void data_to_sketch_1bit(ftr_type o, int dim, int idx);
void write_bit(int offset, int on_off, sketch_type *d);
dist_type dist_L1(ftr_type a, ftr_type b, int dim);
dist_type dist_L1_2(ftr_type a, ftr_type b, int dim, dist_type dist);
dist_type dist_L2(ftr_type a, ftr_type b, int dim);
dist_type dist_L2_2(ftr_type a, ftr_type b, int dim, dist_type dist);
void set_dist_L2_22(ftr_type a);
dist_type dist_L2_22(ftr_type b);
dist_type pdist_L1(ftr_type a, ftr_type b, int dim_from, int dim_to);
dist_type b_dist[MAX_SAMPLESIZE][PJT_DIM];//best
dist_type priority(sketch_type s, query_type *q);
int bit_count(sketch_type a);
void bound(ftr_type q, int b[]);
void bound_with_idx(ftr_type q, bd_type b[]);
void initialize(pivot_type *p);
void get_sample(void);
void resample(int samplesize);
void resample_query(int samplesize);
void init_ind(void);
int check_sample(int x, int s);
int comp(const void *a, const void *b);
int num_sample(int n0, double k1, double k2);
double collision(double *collsize, int idx[], int samplesize);
double collision2(double coll_size, int coll_idx[], int piv);
double correl(unsigned int a[], unsigned int b[], int n);
int comp_sketch(const void *a, const void *b);
int comp_query(const void *a, const void *b);
void print_samplesize(void);
void resample_db(int samplesize);

void print_result(int seed, double db_score, double sc, double time);
int flag = 0;
double Score;
double move = 0,move_s = 0;
int count = 0;
void make_ham_table(void);
void ham_search_16bit(int m);
int comp_sk(const void *a, const void *b);
int comp_qdist(const void *a, const void *b);
int comp_bit(const void *a, const void *b);
int comp_b(const void *a, const void *b);
int comp_bd_type(const void *a, const void *b);
int comp_sk_score(const void *a, const void *b);
int comp_ans(const void *a, const void *b);
void Linf_search_16bit(int m);
int Enum_L1_via_Linf(int k, sketch_type s, bd_type bd_idx[], int sk_num[], sk_score_type seq[]);
void L1_Linf_search_16bit(int m);
void L1_search_16bit_2(int m);
void L1_search_16bit_2_c2_n(int m);
void L1_search_16bit_2_c2_n_new_bucket(int m);
void L1_search_16bit_2_c2_n_new_bucket_para(int m);
void L1_search_16bit_2_kNN(int m);
int countBit16(sketch_type n);
int countBit32(sketch_type n);

static void min_heapify(int i);
int deq(QUE *q);
void enq(QUE *q);

static void min_heapify_p(int i, struct_que *que);
int deq_p(QUE *q, struct_que *que);
void enq_p(QUE *q, struct_que *que);

static void min_heapify_p2(int i, struct_que *que);
int deq_p2(QUE *q, struct_que *que);
void enq_p2(QUE *q, struct_que *que);

void make_empty_que_c(struct_que_c *que);
int new_que_e(struct_que_c *que);
static void min_heapify_c(int i, struct_que_c *que);
int deq_c(int *qe, struct_que_c *que);
void enq_c(int qe, struct_que_c *que);

void make_empty_que_c2(struct_que_c2 *que);
int new_que_e2(struct_que_c2 *que);
static void min_heapify_c2(int i, struct_que_c2 *que);
int deq_c2(QUE_c2 *qe, struct_que_c2 *que);
void deq_c2_del(struct_que_c2 *que);
void enq_c2(QUE_c2 *qe, struct_que_c2 *que);
void enq_c2_and_enq(QUE_c2 *qe, QUE_c2 *qe1, struct_que_c2 *que);
void enq_c2_after_deq(QUE_c2 *qe, struct_que_c2 *que);
void enq_c2_after_deq_and_enq(QUE_c2 *qe, QUE_c2 *qe1, struct_que_c2 *que);

void make_empty_que_c2_n(struct_que_c2_n *que);
int new_que_e2_n(struct_que_c2_n *que);
void min_heapify_c2_n(int i, struct_que_c2_n *que);
int deq_c2_n(QUE_c2 *qe, struct_que_c2_n *que);
void deq_c2_n_del(struct_que_c2_n *que);
void enq_c2_n(QUE_c2 *qe, struct_que_c2_n *que);
void enq_c2_n_after_deq(QUE_c2 *qe, struct_que_c2_n *que);

void make_init_que_c2_s(struct_que_c2_s *que, int num);
void make_empty_que_c2_s(struct_que_c2_s *que);
QUE_c2 *new_que_e2_s(struct_que_c2_s *que);
void min_heapify_c2_s(int i, struct_que_c2_s *que);
int deq_c2_s(QUE_c2 **qe, struct_que_c2_s *que);
void deq_c2_s_del(struct_que_c2_s *que);
void enq_c2_s(struct_que_c2_s *que);
void enq_c2_s_after_deq(struct_que_c2_s *que);

void make_init_que_s(struct_que_s *que, int num);
void make_empty_que_s(struct_que_s *que);
void min_heapify_s(int i, struct_que_s *que);
int deq_s(int *qe, struct_que_s *que);
void deq_s_del(struct_que_s *que);
void enq_s(struct_que_s *que);
void enq_s_after_deq(struct_que_s *que);

void make_empty_que_sk(struct_que_sk *que);
// int new_que_e_sk(struct_que_sk *que);
static void min_heapify_sk(int i, struct_que_sk *que);
int deq_sk(sketch_type *sk, struct_que_sk *que);
void enq_sk(sketch_type sk, struct_que_sk *que);

void search_K_multi(int multi_K[]);
void search_K_multi_c2_n(int multi_K[]);
void search_K_multi_c2_n_new_bucket(int multi_K[]);
//int prepare_bucket_db(int flag);
//void qwrite_bucket_db(char *filename);
//int read_bucket_db(char *filename);
int prepare_bucket_db(int flag, sketch_type *sk, int data_num, int **idx, int **bkt);
void write_bucket_db(char *filename, int data_num, int *idx, int *bkt);
int read_bucket_db(char *filename, int data_num, int **idx, int **bkt);
double precision_resample(int num_K, int resample_size, int mode, double *prec);
double precision_resample_delayed(int num_K, int resample_size, int mode, double *prec);
double precision_resample_cursor(int num_K, int resample_size, int mode, double *prec);
double precision_resample_cursor_2(int num_K, int resample_size, int mode, double *prec);
double precision_resample_cursor_2_n(int num_K, int resample_size, int mode, double *prec);
double precision_resample_cursor_2_s(int num_K, int resample_size, int mode, double *prec);
double precision_resample_s(int num_K, int resample_size, int mode, double *prec);
double precision_resample_cursor_sk(int num_K, int resample_size, int mode, double *prec);
double precision_resample_inf(int num_K, int resample_size, int mode, double *prec);

double precision_resample_by_sequential_filtering(int obj, int resampling_size, int mode, double *prec);


void shuffle(int a[], int n, int m);

#define PARENT(i) ((i)>>1)
#define LEFT(i)   ((i)<<1)
#define RIGHT(i)  (((i)<<1)+1)

#define count_bits(n) (0x0000ffff&((0x00ff00ff&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))\
	+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))\
	+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+\
	+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))+((0xff00ff00&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))\
	+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))\
	+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))>>8)))\
	+((0xffff0000&((0x00ff00ff&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))\
	+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))\
	+((0xff00ff00&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))\
	+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))\
	+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))>>8)))>>16)

#define count_bits2(n) \
	({\
		int x = ((n) & 0x5555) + (((n) >> 1) & 0x5555);\
		x = (x & 0x3333) + ((x >> 2) & 0x3333);\
		x = (x + (x >> 4)) & 0x0F0F;\
		x = (x + (x >> 8)) & 0xFF;\
		x;\
	})

inline int count_bits3( int n){
		n = (n & 0x5555) + ((n >> 1) & 0x5555);
		n = (n & 0x3333) + ((n >> 2) & 0x3333);
		n = (n + (n >> 4)) & 0x0F0F;
		return (n + (n >> 8)) & 0xFF;
}

int mask[32] = {
	1,
	(1 << 1),
	(1 << 2),
	(1 << 3),
	(1 << 4),
	(1 << 5),
	(1 << 6),
	(1 << 7),
	(1 << 8),
	(1 << 9),
	(1 << 10),
	(1 << 11),
	(1 << 12),
	(1 << 13),
	(1 << 14),
	(1 << 15),
	(1 << 16),
	(1 << 17),
	(1 << 18),
	(1 << 19),
	(1 << 20),
	(1 << 21),
	(1 << 22),
	(1 << 23),
	(1 << 24),
	(1 << 25),
	(1 << 26),
	(1 << 27),
	(1 << 28),
	(1 << 29),
	(1 << 30),
	(1 << 31)
};

#ifdef USE_MASK_PATTERN
#define MASK(i) (mask[i])
#else
#define MASK(i) (1 << (i))
#endif

static void min_heapify(int i)
{
    int l, r;
    int smallest;

    l = LEFT(i);
	r = RIGHT(i);
    if (l < qsize && que[l].key < que[i].key) smallest = l; else smallest = i;
    if (r < qsize && que[r].key < que[smallest].key) smallest = r;
    if (smallest != i) {
        QUE t = que[i]; que[i] = que[smallest]; que[smallest] = t;
        min_heapify(smallest);
    }
}

int deq(QUE *q)
{
    if (qsize == 0) return 0;
    memcpy(q, &que[0], sizeof(QUE));
    que[0] = que[--qsize];
    min_heapify(0);
    return 1;
}

void enq(QUE *q)
{
    int i, ii;

    i = qsize++;
    memcpy(&que[i], q, sizeof(QUE));
	while (i > 0 && que[ii = PARENT(i)].key > que[i].key) {
        QUE t = que[i];
		que[i] = que[ii]; 
		que[ii] = t;
        i = ii;
    }
}

// enq and deq for parallelize
static void min_heapify_p(int i, struct_que *que)
{
    int l, r;
    int smallest;

    l = LEFT(i);
	r = RIGHT(i);
    if (l < que->qsize && que->element[l].key < que->element[i].key) smallest = l; else smallest = i;
    if (r < que->qsize && que->element[r].key < que->element[smallest].key) smallest = r;
    if (smallest != i) {
        QUE t = que->element[i]; que->element[i] = que->element[smallest]; que->element[smallest] = t;
        min_heapify_p(smallest, que);
    }
}

int deq_p(QUE *q, struct_que *que)
{
    if (que->qsize == 0) return 0;
    memcpy(q, &(que->element[0]), sizeof(QUE));
	que->element[0] = que->element[--(que->qsize)];
    min_heapify_p(0, que);
    return 1;
}

void enq_p(QUE *q, struct_que *que)
{
    int i, ii;

	i = (que->qsize)++;
    memcpy(&(que->element[i]), q, sizeof(QUE));
	while (i > 0 && que->element[ii = PARENT(i)].key > que->element[i].key) {
        QUE t = que->element[i];
		que->element[i] = que->element[ii]; 
		que->element[ii] = t;
        i = ii;
    }
}

#define LESS_P(a, b) (que->element[(a)].key < que->element[(b)].key || (que->element[(a)].key == que->element[(b)].key && que->element[(a)].sk < que->element[(b)].sk))
#define LARG_P(a, b) (que->element[(a)].key > que->element[(b)].key || (que->element[(a)].key == que->element[(b)].key && que->element[(a)].sk > que->element[(b)].sk))

// enq and deq for parallelize ver Imamura
static void min_heapify_p2(int i, struct_que *que)
{
    QUE t = que->element[i];

	int temp = QSIZE;
	que->element[temp] = t;

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (que->element[l].key < t.key) {
            int r = RIGHT(i);
            if (r < que->qsize && LESS_P(r, l)) {
                que->element[i] = que->element[r];
                i = r;
            }
            else {
                que->element[i] = que->element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && LESS_P(r, temp))) {
                break;
            }
            que->element[i] = que->element[r];
            i = r;
        }
    }
    que->element[i] = t;
}

int deq_p2(QUE *q, struct_que *que)
{
    if (que->qsize == 0) return 0;
    *q = que->element[0];
    que->element[0] = que->element[--(que->qsize)];
    min_heapify_p2(0, que);
    return 1;
}

void enq_p2(QUE *q, struct_que *que)
{
    int i, ii;

    i = (que->qsize)++;
    que->element[i] = *q;
    while (i > 0 && LARG_P(PARENT(i), i)) {
    	ii = PARENT(i);
        QUE t = que->element[i];
        que->element[i] = que->element[ii]; 
        que->element[ii] = t;
        i = ii;
    }
}

// enq and deq for parallelize using cursors by Takeshi
void make_empty_que_c(struct_que_c *que)
{
	que->qsize = 0;
	que->used = 0;
}

int new_que_e(struct_que_c *que)
{
	int i = que->used;
	que->used++;
	return i;
}

#define LESS_C(a, b) (que_e[element[(a)]].key < que_e[element[(b)]].key || (que_e[element[(a)]].key == que_e[element[(b)]].key && que_e[element[(a)]].sk < que_e[element[(b)]].sk))

static void min_heapify_c(int i, struct_que_c *que)
{
	que_elm *que_e = que->que_e;
	int *element = que->element;
    int t = element[i];
	
	int temp = QSIZE;
	element[temp] = t;

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (que_e[element[l]].key < que_e[t].key) {
            int r = RIGHT(i);
            if (r < que->qsize && LESS_C(r, l)) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && LESS_C(r, temp))) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_c(int *qe, struct_que_c *que)
{
    if (que->qsize == 0) return 0;
    *qe = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c(0, que);
    return 1;
}

void enq_c(int qe, struct_que_c *que)
{
	que_elm *que_e = que->que_e;
	int *element = que->element;
    int i, ii;
	int t;

    i = (que->qsize)++;
    element[i] = qe;
    while (i > 0 && LESS_C(i, PARENT(i))) {
    	ii = PARENT(i);
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

// enq and deq for parallelize using cursors by Yasunobu
void make_empty_que_c2(struct_que_c2 *que)
{
	que->qsize = 0;
	que->detail_size = 0;
}

int new_que_e2(struct_que_c2 *que)
{
	int i = que->detail_size;
	que->detail_size++;
	return i;
}

static void min_heapify_c2(int i, struct_que_c2 *que)
{
	QUE_c2 *element = que->element;
    QUE_c2 t = element[i];

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (element[l].key < t.key) {
            int r = RIGHT(i);
            if (r < que->qsize && element[r].key < element[l].key) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && element[r].key < t.key)) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_c2(QUE_c2 *qe, struct_que_c2 *que)
{
    if (que->qsize == 0) return 0;
	#ifdef USE_TOP
	*qe = que->element[0]; // Top 要素を返すだけで，Queue 本体では，Top要素を取り除かない
	return 1;
	#else
    *qe = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c2(0, que);
    return 1;
	#endif
}

void deq_c2_del(struct_que_c2 *que) // Deq で取り除かなかったものを取り除く 
{
    if (que->qsize == 0) return;
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c2(0, que);
}

void enq_c2(QUE_c2 *qe, struct_que_c2 *que)
{
	QUE_c2 *element = que->element;
    int i, ii;
	QUE_c2 t;

    i = (que->qsize)++;
    element[i] = *qe;
    while (i > 0 && element[ii = PARENT(i)].key > element[i].key) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void enq_c2_and_enq(QUE_c2 *qe, QUE_c2 *qe1, struct_que_c2 *que)
{
	QUE_c2 *element = que->element;
    int i, ii;
	QUE_c2 t;

    i = (que->qsize)++;
    element[i] = *qe;
    while (i > 0 && element[ii = PARENT(i)].key > element[i].key) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }

    i = (que->qsize)++;
    element[i] = *qe1;
    while (i > 0 && element[ii = PARENT(i)].key > element[i].key) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void enq_c2_after_deq(QUE_c2 *qe, struct_que_c2 *que)
{
    que->element[0] = *qe;
    min_heapify_c2(0, que);
}

// enq and deq for parallelize using cursors by Yasunobu (modified for Naoya's version)
#define CURSOR_STACK_SIZE 40000
void make_empty_que_c2_n(struct_que_c2_n *que)
{
	que->qsize = 0;
	que->detail_size = 0;
	#ifdef CURSOR_STACK
	for(int i = 0; i < CURSOR_STACK_SIZE; i++) que->cursor_stack[i] = CURSOR_STACK_SIZE - 1 - i;
	que->stack_top =    CURSOR_STACK_SIZE - 1;
	#endif
}

int new_que_e2_n(struct_que_c2_n *que)
{
	#ifdef CURSOR_STACK
	int i = que->cursor_stack[que->stack_top--];
	#else
	int i = que->detail_size;
	que->detail_size++;
	#endif
	return i;
}

void min_heapify_c2_n(int i, struct_que_c2_n *que)
{
	QUE_c2 *element = que->element;
    QUE_c2 t = element[i];

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (element[l].key < t.key) {
            int r = RIGHT(i);
            if (r < que->qsize && element[r].key < element[l].key) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && element[r].key < t.key)) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_c2_n(QUE_c2 *qe, struct_que_c2_n *que)
{
    if (que->qsize == 0) return 0;
	#ifdef USE_TOP
    *qe = que->element[0];
	return 1;
	#else
    *qe = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c2_n(0, que);
    return 1;
	#endif
}

void deq_c2_n_del(struct_que_c2_n *que) // Deq で取り除かなかったものを取り除く 
{
    if (que->qsize == 0) return;
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c2_n(0, que);
}

void enq_c2_n(QUE_c2 *qe, struct_que_c2_n *que)
{
	QUE_c2 *element = que->element;
    int i, ii;
	QUE_c2 t;

    i = (que->qsize)++;
    element[i] = *qe;
    while (i > 0 && element[ii = PARENT(i)].key > element[i].key) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void enq_c2_n_after_deq(QUE_c2 *qe, struct_que_c2_n *que)
{
    que->element[0] = *qe;
    min_heapify_c2_n(0, que);
}

// enq and deq for parallelize using cursors by Yasunobu (modified for Naoya's version)
// using Heap (shino)
void make_init_que_c2_s(struct_que_c2_s *que, int num)
{
	que->qsize = 0;
	for(int i = 0; i < num; i++) {
		que->element[i].cursor = i;
	}
}

void make_empty_que_c2_s(struct_que_c2_s *que)
{
	que->qsize = 0;
}

QUE_c2 *new_que_e2_s(struct_que_c2_s *que)
{
	return &(que->element[que->qsize]);
}

void min_heapify_c2_s(int i, struct_que_c2_s *que)
{
	QUE_c2 *element = que->element;
    QUE_c2 t = element[i];

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (element[l].key < t.key) {
            int r = RIGHT(i);
            if (r < que->qsize && element[r].key < element[l].key) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && element[r].key < t.key)) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_c2_s(QUE_c2 **qe, struct_que_c2_s *que)
{
    if (que->qsize == 0) return 0;
	*qe = &(que->element[0]);
	return 1;
}

void deq_c2_s_del(struct_que_c2_s *que) // Deq で取り除かなかったものを取り除く 
{
	QUE_c2 t;
    if (que->qsize == 0) return;
	t = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
    que->element[que->qsize - 1] = t;
	--(que->qsize);
    min_heapify_c2_s(0, que);
}

void enq_c2_s(struct_que_c2_s *que)
{
	QUE_c2 *element = que->element;
    int i, ii;
	QUE_c2 t;

    i = (que->qsize)++;
    while (i > 0 && element[ii = PARENT(i)].key > element[i].key) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void enq_c2_s_after_deq(struct_que_c2_s *que)
{
    min_heapify_c2_s(0, que);
}

// using Heap as Stack (shino)
void make_init_que_s(struct_que_s *que, int num)
{
	que->qsize = 0;
	for(int i = 0; i < num; i++) {
		que->element[i] = i;
	}
}

void make_empty_que_s(struct_que_s *que)
{
	que->qsize = 0;
}

#define KEY(i) (que->details[i].key)
#define SK(i) (que->details[i].sk)
#define PT(i) (que->details[i].pt)

void min_heapify_s(int i, struct_que_s *que)
{
	int *element = que->element;
    int t = element[i];

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (KEY(element[l]) < KEY(t)) {
            int r = RIGHT(i);
            if (r < que->qsize && KEY(element[r]) < KEY(element[l])) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && KEY(element[r]) < KEY(t))) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_s(int *qe, struct_que_s *que)
{
    if (que->qsize == 0) return 0;
	*qe = que->element[0];
	return 1;
}

void deq_s_del(struct_que_s *que) // Deq で取り除かなかったものを取り除く 
{
	int t;
    if (que->qsize == 0) return;
	t = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
    que->element[que->qsize - 1] = t;
	--(que->qsize);
    min_heapify_s(0, que);
}

void enq_s(struct_que_s *que)
{
	int *element = que->element;
    int i, ii;
	int t;

    i = (que->qsize)++;
    while (i > 0 && KEY(element[ii = PARENT(i)]) > KEY(element[i])) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void enq_s_after_deq(struct_que_s *que)
{
    min_heapify_s(0, que);
}

// enq and deq for parallelize using cursors by Takeshi
// using sketches as cursors
void make_empty_que_sk(struct_que_sk *que)
{
	que->qsize = 0;
}

static void min_heapify_sk(int i, struct_que_sk *que)
{
	dist_type *key = que->key;
	int *element = que->element;
    sketch_type t = element[i];
	
	int temp = QSIZE;
	element[temp] = t;

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (key[element[l]] < key[t]) {
            int r = RIGHT(i);
            if (r < que->qsize && key[element[r]] < key[element[l]]) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && key[element[r]] < key[t])) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_sk(sketch_type *sk, struct_que_sk *que)
{
    if (que->qsize == 0) return 0;
    *sk = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_sk(0, que);
    return 1;
}

void enq_sk(sketch_type sk, struct_que_sk *que)
{
//	que_elm_sk *que_e = que->que_e;
	dist_type *key = que->key;
	int *element = que->element;
    int i, ii;
	sketch_type t;

    i = (que->qsize)++;
    element[i] = sk;
    while (i > 0 && key[element[ii = PARENT(i)]] > key[element[i]]) {
    	ii = PARENT(i);
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void print_result(int seed, double db_score, double sc, double time){
	
	FILE *outputfile;
	outputfile = fopen("allresult.txt", "a+");  

	if(outputfile == NULL){
		printf("cannot open\n");
		exit(1);                      
	}
	fprintf(outputfile, "------%d------\n",seed);
	fprintf(outputfile, "%.0f\n%.0f\n%.0f\n%.2f\n",Score, db_score, sc, time); 
	fclose(outputfile);
}

void init_ind(void){
	int i;
	for(i = 0; i < MAX_SAMPLESIZE; i++){
		ind_re[i] = i;
	}

}

int comp(const void *a, const void *b) {
	
	if(*((unsigned int *) a) < *((unsigned int *) b))
		return -1;
	else if(*((unsigned int *) a) == *((unsigned int *) b))
		return 0;
	else
		return 1;
}

int comp_sk(const void *a, const void *b) {
	
	
	if(sketch[*((int *) a)] < sketch[*((int *) b)])
		return -1;
	else if(sketch[*((int *) a)] == sketch[*((int *) b)])
		return 0;
	else
		return 1;
}

int comp_qdist(const void *a, const void *b) {
	if(((search_type *) a) -> dist < (((search_type *) b) -> dist))
		return -1;
	else if(((search_type *) a) -> dist == (((search_type *) b) -> dist))
		return 0;
	else
		return 1;
}

int comp_ans(const void *a, const void *b) {
	if(((answer *) a) -> dist < (((answer *) b) -> dist))
		return -1;
	else if(((answer *) a) -> dist == (((answer *) b) -> dist))
		return 0;
	else
		return 1;
}

int comp_bit(const void *a, const void *b) {
	if(bitcnt_tbl[*((int *) a)] < bitcnt_tbl[*((int *) b)])
		return -1;
	else if(bitcnt_tbl[*((int *) a)] == bitcnt_tbl[*((int *) b)])
		return 0;
	else
		return 1;
}

int comp_b(const void *a, const void *b) {
	if(bd[*((int *) a)] < bd[*((int *) b)])
		return -1;
	else if(bd[*((int *) a)] == bd[*((int *) b)])
		return 0;
	else
		return 1;
}

int comp_bd_type(const void *a, const void *b) {
	if(((bd_type *) a)->bd < ((bd_type *) b)->bd)
		return -1;
	else if(((bd_type *) a)->bd == ((bd_type *) b)->bd)
		return 0;
	else
		return 1;
}

int comp_sk_score(const void *a, const void *b) {
	if(((sk_score_type *) a)->score < ((sk_score_type *) b)->score)
		return -1;
	else if(((sk_score_type *) a)->score == ((sk_score_type *) b)->score)
		return 0;
	else
		return 1;
}

void get_sample(void) {
	int i;
	
	for(i = 0; i < MAX_SAMPLESIZE;){
		
		sample[i] = (ftr_type)calloc(db_header->data_dim, db_header->element_size);
		
		int x = random() % db_header->data_num;
		int j;

		if(check_sample(x, i) == 1) {
			for(j = 0; j < db_header->data_dim; j++) {
				sample[i][j] = database[x][j];
			}
			sample_number[i] = x; // sample_numberに登録
			i++;
		}
	 }
#if SKETCH_OPTIMIZATION == 1
// 2i番目と2i+1番目のsample間の実距離を求めておく
// 射影距離も0に初期化しておく
	for(i = 0; i < MAX_SAMPLESIZE / 2; i++){
		od[i] = DISTANCE(sample[2 * i], sample[2 * i + 1], db_header->data_dim);
		b_d1[i] = 0;
	}
#endif
	//fclose(fp);
}

void resample_db(int samplesize) {
	int i;
	int x;
	ftr_type temp;

	for(i = 0; i < samplesize; i++){
			x = random() % (db_header->data_num-i) + i;
			temp = database_temp[i];
			database_temp[i] = database_temp[x];
			database_temp[x] = temp;
	}
}

void resample(int samplesize) {
	int i,x;
	ftr_type temp;
	sketch_type sk_temp;
	
		for(i = 0; i < samplesize; i++){
			x = random() % (MAX_SAMPLESIZE-i) + i;
			temp = sample[i];
			sample[i] = sample[x];
			sample[x] = temp;
			sk_temp = sample_sketch[i];
			sample_sketch[i] = sample_sketch[x];
			sample_sketch[x] = sk_temp;
		}
}

void resample_query(int samplesize) {
	int i,x;
	query_type temp;
	
	for(i = 0; i < samplesize; i++){
		x = random() % (query_num - i) + i;
		temp = query[i];
		query[i] = query[x];
		query[x] = temp;
	}
//	qsort(query, samplesize, sizeof(query_type), comp_query);

}

int check_sample(int x, int s){
	int i;

	for(i = 0; i < s; i++){
		if(sample_number[i] == x){
			return 0;
		}
	}
	return 1;
}

dist_type compute_rad(ftr_type piv, int dim, int samplesize, int *count_tbl) {
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

void initialize(pivot_type *p){
	int dim, axi;
	for(dim = 0; dim < PJT_DIM; dim++) {
		p->piv0[dim] = (ftr_type)malloc(db_header->element_size * db_header->data_dim);
#if SKETCH_PARTITION == 0
		p->piv1[dim] = (ftr_type)malloc(db_header->element_size * db_header->data_dim);
#endif
		for(axi = 0; axi < db_header->data_dim; axi++) {
			p->piv0[dim][axi] = 0;
#if SKETCH_PARTITION == 0
			p->piv1[dim][axi] = 0;
#endif
			p->rad[dim] = 0;
		}
	}
	p->partition = 0;
#ifndef NO_SAMPLE_SKETCH
	sample_sketch = (sketch_type *)calloc(MAX_SAMPLESIZE , sizeof(sketch_type));
#endif
}

void print_bin(unsigned int s)
{
	int i;
	for(i = 0; i < 32; i++)
		printf("%1d",(s >> (31 - i)) & 1);
}

int msb_pos(unsigned int x) {
  int pos = -1;
 
  if (x != 0) {
    __asm__("bsrl %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

int lsb_pos(unsigned int x) {
  int pos = -1;
 
  if (x != 0) {
    __asm__("bsfl %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

int msb_pos_long(unsigned long x) {
  long pos = -1;
 
  if (x != 0) {
    __asm__("bsrq %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

int lsb_pos_long(unsigned long x) {
  long pos = -1;
 
  if (x != 0) {
    __asm__("bsfq %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

#ifndef PART_TERM
#define PART_TERM 4
#endif

dist_type dist_L1(ftr_type a, ftr_type b, int dim)
// 距離関数（L1）
{
	int j;
	dist_type s = 0;
	
	for(j = 0; j < dim; j++)  {
		s += abs((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type dist_L1_2(ftr_type a, ftr_type b, int dim, dist_type dist)
// 距離関数（L1）
{
	int j, k;
	dist_type s = 0;
	
	for(k = 0; k < PART_TERM; k++) {
		for(j = k * dim / PART_TERM; j < (k + 1) * dim / PART_TERM; j++) {
			s += abs((int)a[j] - (int)b[j]);
		}
		if(s >= dist) break;
	}
	return s;
}

dist_type dist_L2(ftr_type a, ftr_type b, int dim)
// 距離関数（L2）
{
	int j;
	dist_type s = 0;
	
	for(j = 0; j < dim; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type dist_L2_2(ftr_type a, ftr_type b, int dim, dist_type dist)
// 距離関数（L2）
{
	int j, k;
	dist_type s = 0;
	
	for(k = 0; k < PART_TERM; k++) {
		for(j = k * dim / PART_TERM; j < (k + 1) * dim / PART_TERM; j++) {
			s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
		}
		if(s >= dist) return s;
	}
	return s;
}

dist_type dist_L2_21(ftr_type a, ftr_type b)
// 距離関数（L2）
{
	int j;
	dist_type s = 0;
	
	for(j = 0; j < FTR_DIM; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return s;
}

int point_a[FTR_DIM];

void set_dist_L2_22(ftr_type a)
{
	int j;
	for(j = 0; j < FTR_DIM; j++) 
		point_a[j] = a[j];
}

dist_type dist_L2_22(ftr_type b)
// 距離関数（L2）
{
	int j;
	dist_type s = 0;
	for(j = 0; j < FTR_DIM; j++)  {
		s += (point_a[j] - (int)b[j]) * (point_a[j] - (int)b[j]);
	}
	return s;
}

dist_type pdist_L1(ftr_type a, ftr_type b, int dim_from, int dim_to)
// 座標分割法のための部分空間の距離関数（L1）
// ftr_type a, b: 特徴データ，dim: 特徴データの次元数，
// dim_from, dim_to: 部分空間の軸の最初と終わり
{
	int j;
	dist_type s = 0;
	
	for(j = dim_from; j <= dim_to; j++)  {
		s += abs((int)a[j] - (int)b[j]);
	}
	return s;
}

#ifdef EVAL_BY_SEQUENTIAL_FILTERING
dist_type priority(sketch_type s, query_type *q)
{
	sketch_type d = s ^ q->sketch;
	return q->tbl[0][d & 0xff] + q->tbl[1][(d >> 8) & 0xff] + q->tbl[2][(d >> 16) & 0xff] + q->tbl[3][d >> 24];
}
#endif

void write_bit(int offset, int on_off, sketch_type *d)
{
	if(on_off == 1) {
		*d = *d | (1 << offset);
	} else {
		*d = *d & ~(1 << offset);
	}
}

void median(void)
{
	int cnt[256];
	int i, j, s, t;
	
	med = (ftr_type)malloc(db_header->element_size * db_header->data_dim);

	for(i = 0; i < 256 /* BIT */; i++) {
		cnt[i] = 0;
	}
	for(j = 0; j < db_header->data_dim; j++) {
		for(i = t = 0; i < db_header->data_num; i += 500) {
			cnt[database[i][j]]++; t++;
		}
	}
	for(j = 0; j < db_header->data_dim; j++) {
		for(i = 0; i < 256 /* BIT */; i++) {
			cnt[i] = 0;
		}
		
		for(i = t = 0; i < db_header->data_num; i += 500) {
			cnt[database[i][j]]++; t++;
		}
		for(i = s = 0; (i < 256 /* BIT */) && (s < t / 2); i++) {
			s += cnt[i];
		}
		med[j] = i;
	}
}

int num_sample(int n0, double k1, double k2){
	double remain = 1.0 - progress;
	double n = n0*(1/pow((pow(remain,k1)*pow(2, -progress*k2)),2));
	if(n > MAX_SAMPLESIZE)
		return MAX_SAMPLESIZE;
	else
		return (int)n;
}

int difference_QBP(ftr_type c1, ftr_type c2) {
	int i, d = 0;
	for(i = 0; i < FTR_DIM; i++) {
		d += (c1[i] != c2[i]);
	}
	return d;
}

void make_pivot_for_QBP(int c, ftr_type center, int dim) {
	int i = 0;
	int dim_from = 0, dim_to = 0;
#if SKETCH_PARTITION == 3
	dim_from = dim * 4;
	dim_to = (dim + 1) * 4 - 1;
#endif
//	c = random() % db_header->data_num; // データベースから乱択

	if(pivot.partition == 1)   { // BP
		for(i = 0; i < db_header->data_dim; i++) {
			pivot.piv0[dim][i] = database[c][i];
		}
		memcpy(center, database[c], db_header->element_size * db_header->data_dim);
	} else if(pivot.partition == 2)   { // QBP
		for(i = 0; i < db_header->data_dim; i++) {
			center[i] = (database[c][i] < med[i]) ? 0 : 255;
		}
	} else if(pivot.partition == 3)   { // PQBP
		for(i = dim_from; i <= dim_to; i++) {
			center[i] = (database[c][i] < med[i]) ? 0 : 255;
		}
	}
}

void select_pivot_BPorQBP(void) {
	int i = 0, dim, t;
#if SKETCH_PARTITION == 3
	int dim_from = 0, dim_to = 0;
#endif
#if SKETCH_OPTIMIZATION == 0
	double scmin = DBL_MAX;
#else
	double scmax = -DBL_MAX;
#endif
	double coll_size;
	int coll_idx[MAX_SAMPLESIZE];
	int count_tbl[MAX_SAMPLESIZE];
	dist_type r = 0;
	ftr_type max_piv = (ftr_type)malloc(db_header->element_size * db_header->data_dim);

	for(dim = 0; dim < PJT_DIM; dim++){	

		//scmin = DBL_MAX;
#if SKETCH_OPTIMIZATION == 0
		scmin = collision(&coll_size, coll_idx, MAX_SAMPLESIZE);
//		fprintf(stderr, "scmin = %8lf\n", scmin);

#elif SKETCH_OPTIMIZATION == 1
		scmin = -correl();
#endif

#if SKETCH_PARTITION == 3
		dim_from = dim * 4;
		dim_to = (dim + 1) * 4 - 1;
#endif

		for(t = 0; t < NUM_TRIAL1; t++){	// NUM_TRIAL=試行回数 
			make_pivot_for_QBP(random() % db_header->data_num, pivot.piv0[dim], dim);
			pivot.rad[dim] = compute_rad(pivot.piv0[dim], db_header->data_dim, MAX_SAMPLESIZE / 10, count_tbl);  // 半径計算 (BP と同じ計算方法の場合)
			double sc = collision2(coll_size, coll_idx, dim);
			if(sc < scmin){
				for(i = 0; i < db_header->data_dim; i++)
					max_piv[i] = pivot.piv0[dim][i];
				scmin = sc;
				r = pivot.rad[dim];
			}
		}//NUM_TRIAL終了
		
		for(i = 0; i < db_header->data_dim; i++)
			pivot.piv0[dim][i] = max_piv[i]; //best pivot
		pivot.rad[dim] = r;                  //best radian
		make_sample_sketch(dim, MAX_SAMPLESIZE);    // Bestなpivotでsketchを作り直す（Trial中にsketchを書き換えているため）
		fprintf(stderr,"DIM%3d fin ... = %12.0f\n", dim+1, scmin);
	}
	fprintf(stderr,"DIM%3d fin ... = %12.0f, %12.4e\n", dim+1, scmin, scmin / (MAX_SAMPLESIZE * (MAX_SAMPLESIZE - 1) / 2));
	fprintf(stdout,"DIM%3d fin ... = %12.0f, %12.4e\n", dim+1, scmin, scmin / (MAX_SAMPLESIZE * (MAX_SAMPLESIZE - 1) / 2));
	
	free(max_piv);
}

void random_QBP(void) {
	int i, dim;

	if(pivot.partition != 2) {
		fprintf(stderr, "Sorry, random_QBP is only for QBP\n");
		exit(1);
	}
	for(dim = 0; dim < PJT_DIM; dim++){	
		pivot.rad[dim] = 0;
		for(i = 0; i < db_header->data_dim; i++) {
			pivot.piv0[dim][i] = (random() % 2) ? 0 : 255;
			pivot.rad[dim] += abs(pivot.piv0[dim][i] - med[i]); 
		}
	}
}

#ifndef EVAL_PRECISION
#define EVAL_PRECISION precision_resample
#endif

int debug_flag = 0;

void optimize_pivot_by_precision_AIR_flip(int obj, int mode, int f_num) {
	// obj: 目的関数の指定（精度）   1 <= obj <= 69 → K (%) = obj / 10 で検索したときの精度）
	//                    （候補数）90 <= obj < 100 → 精度が obj (%) になる候補数 K(%)
	// ただし，現状では上の精度でのみデバッグ・調整中
	// mode: 評価の集計方法（0: 正解数の割合，1: 正解の重み付き合計）
	// f_num:  1回の試行で同時にフリップ（0 と 255 を反転）する次元数
	// 【注意】とりあえず，座標分割法には未対応なので，SKETCH_PARTITION = 3 では使用しないこと
	// 質問集合に対して，再サンプリングの最終サイズを半分で止めるようにする．#define USE_HALF 
	// USE_HALF は精度を目的関数とするときのみ，つまり，1 <= obj <= 19 にのみ対応
	int i = 0;
	int dim, axi[FTR_DIM];
	int multi_K[1000], accept;
	double eval, eval_best, prec;
	int n0 = N0;
	int reuse = REUSE, reuse_c = 0;
	int resampling_size = 10;
	double Tr;
#ifdef VAR_FLIP
	#ifndef REPLACE_WHOLE
	int f_max = f_num; // フリップの最大数
	#endif
#endif
	int not_improved = 0, final_check = 0;
#ifndef CONVERGENCE_CHECK
#define CONVERGENCE_CHECK 0
#endif
//	int start, end;
	int max_query_num = query_num;
//	sketch_type *temp_sketch;
	int num_accept = 0;
#ifdef USE_HALF
	max_query_num /= 2;
#endif
#ifdef REPLACE_WHOLE
	#ifndef REPLACE_TRIAL
		#define REPLACE_TRIAL 1
	#endif
	int j;
	int candidate[REPLACE_TRIAL];
	static ftr_type piv_candidate[REPLACE_TRIAL] = {NULL};
	static ftr_type piv_temp = NULL;
	dist_type rad_temp;
	int r_count;
	int piv_diff[REPLACE_TRIAL], piv_diff_min, piv_diff_min_n = 0;
	int count_tbl[MAX_SAMPLESIZE / 10];
#else
	int count_tbl[MAX_SAMPLESIZE / 10];
#endif

//	max_query_num = query_num = 1000;

	if(SKETCH_PARTITION == 3) {
		fprintf(stderr, "Sorry!  optimize_pivot_by_precision_AIR_flip is not applicable to PQBP\n");
		exit(1);
	}
	
	ans = (answer *)malloc((query_num * 2) * sizeof(answer));

	make_sketch();
	if(y_sketch == NULL) {
		y_sketch = (sketch_type *)malloc(sizeof(sketch_type) * db_header->data_num);
	}
//	make_query_sketch();
	make_query_sketch_for_resample(max_query_num);
	init_tables();
	make_ham_table();
	for(i = 0; i < FTR_DIM; i++) {
		axi[i] = i;
	}
	if(obj > 0 && obj < 90) {
		eval_best = EVAL_PRECISION(obj, max_query_num, mode, &prec);
		fprintf(stderr, "(1) Precision for K = %.1lf%% = %6.4lf%%, reuse = %d\n", (double)obj / 10, eval_best, reuse);
	} else {
		search_K_multi(multi_K);
		eval_best = multi_K[obj * 10];
	}

	for(i = 0; i < 1000; i++) {
		resample_query(max_query_num);
	}
	
	for(i = reuse_c = 0; i < NUM_TRIAL2 + NUM_LS + CONVERGENCE_CHECK; i++, reuse_c = (reuse_c + 1) < reuse ? reuse_c + 1 : 0) {
		if(i == NUM_TRIAL2) { // DO Local Search using whole queries
			resampling_size = max_query_num;
			fprintf(stderr, "Goto local search: trials = %d\n", NUM_LS);
			printf("Goto local search: trials = %d\n", NUM_LS);
			make_query_sketch_for_resample(resampling_size);
			eval_best = EVAL_PRECISION(obj, resampling_size, mode, &prec);
			fprintf(stderr, "(2) Precision for K = %.1lf%% = %6.4lf%%, size = %5d, i = %6d\n", (double)obj / 10, eval_best, resampling_size, i);
		}
		if(i == NUM_TRIAL2 + NUM_LS) {
			fprintf(stderr, "Goto convergence check: trials = %d\n", CONVERGENCE_CHECK);
			printf("Goto convergence check: trials = %d\n", CONVERGENCE_CHECK);
			final_check = 1;
		} else if(i < NUM_TRIAL2 && reuse_c == 0) {
			int Ti = i + REUSE;
			if(Ti > NUM_TRIAL2) Ti = NUM_TRIAL2 - 1;
			Tr = 1 - (double)Ti / (NUM_TRIAL2 == 0 ? 1 : NUM_TRIAL2); // Tr = remain; 1 - Tr = progress;
			if(T_POWER >= 1.0) {
				Tr = pow(Tr, T_POWER);
			} else {
				Tr = Tr * pow(T_POWER, 1 - Tr); // remain^aplha * beta^progress; alpha = 1, beta = T_POWER;
			}
			resampling_size = (double)max_query_num / ((double)(max_query_num - n0) / n0 * Tr * Tr + 1);
			reuse_c = 0;
			resample_query(resampling_size);
			make_query_sketch_for_resample(resampling_size);
			eval_best = EVAL_PRECISION(obj, resampling_size, mode, &prec);
			fprintf(stderr, "(3) Precision for K = %.1lf%% = %6.4lf%%, size = %5d, i = %6d\r", (double)obj / 10, eval_best, resampling_size, i);
			
		}
		// ピボット（dim = 0, ... , PJT_DIM)と次元 (axi = 0, ... , FTR_DIM) を乱択
		// pivot.piv0[dim][axi] の値を反転（フリップ）（0 <==> 255)
		// 半径を更新
		if(final_check) {
			if(final_check == 1) {
				dim = 0;
				axi[0] = 0;
				f_num = 1;
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
			fprintf(stderr, "Convergence check (i = %7d): dim = %2d, axi = %2d, not_improved = %4d\r", i - NUM_TRIAL2 - NUM_LS, dim, axi[0], not_improved);
		} else {
			dim = random() % PJT_DIM;
#ifdef REPLACE_WHOLE
// ピボットを一つ(dim次元のもの)を全部取り換える．
// 現在のピボットを退避する．
			if(piv_candidate[0] == NULL) {
				for(j = 0; j < REPLACE_TRIAL; j++) { 
					piv_candidate[j] = (ftr_type)malloc(db_header->element_size * db_header->data_dim);
				}
			}
			for(r_count = 0; r_count < REPLACE_TRIAL; r_count++) {
				candidate[r_count] = random() % db_header->data_num;
			}
			#pragma omp parallel for
			for(r_count = 0; r_count < REPLACE_TRIAL; r_count++) {
				make_pivot_for_QBP(candidate[r_count], piv_candidate[r_count], dim);
				piv_diff[r_count] = difference_QBP(piv_candidate[r_count], pivot.piv0[dim]);
			}
			piv_diff_min = 1000;
			for(r_count = 0; r_count < REPLACE_TRIAL; r_count++) {
				if(piv_diff[r_count] < piv_diff_min) {
					piv_diff_min = piv_diff[r_count];
					piv_diff_min_n = r_count;
				}
			}
			if(piv_temp == NULL) {
				piv_temp = (ftr_type)malloc(db_header->element_size * db_header->data_dim);
			}
			memcpy(piv_temp, pivot.piv0[dim], db_header->element_size * db_header->data_dim);
			rad_temp = pivot.rad[dim];
			memcpy(pivot.piv0[dim], piv_candidate[piv_diff_min_n], db_header->element_size * db_header->data_dim);
			pivot.rad[dim] = compute_rad(piv_candidate[piv_diff_min_n], db_header->data_dim, MAX_SAMPLESIZE / 10, count_tbl);  // 半径計算 (BP と同じ計算方法の場合)
		}
#else // REPLACE_WHOLE
	#ifdef VAR_FLIP
			f_num = (int)(0.9 + sqrt((double)(NUM_TRIAL2 + NUM_LS - i) / (NUM_TRIAL2 + NUM_LS)) * f_max);
			if(f_num <= 0) f_num = 1;
	#endif
			shuffle(axi, FTR_DIM, f_num);

		}
		for(int j = 0; j < f_num; j++) {
			pivot.piv0[dim][axi[j]] = 255 - pivot.piv0[dim][axi[j]];	// 0 <=> 255 のフリップ
		}
		pivot.rad[dim] = compute_rad(pivot.piv0[dim], db_header->data_dim, MAX_SAMPLESIZE / 10, count_tbl);  // 半径計算 (BP と同じ計算方法の場合)
#endif
		make_y_sketch(dim);
		remake_query_sketch_for_resample(dim, resampling_size);
		accept = 0;
		if(obj > 0 && obj < 90) {
			temp_sketch = sketch; sketch = y_sketch; y_sketch = temp_sketch; // 精度評価のために sketch と y_sketch を交換
			eval = EVAL_PRECISION(obj, resampling_size, mode, &prec);
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
			} else if(final_check) {
				not_improved++;
			}
			if(i % 100 == 0) {
				printf("(4) Precision for K = %.1lf%% = %10.8lf%%, size = %5d, i = %6d\n", (double)obj / 10, eval_best, resampling_size, i);
				#ifdef INTERMEDIATE
					eval = EVAL_PRECISION(obj, max_query_num, mode, &prec);
					printf("%6d, %6.2lf, INTERMEDIATE\n", i, eval);
				#endif
			}
		} else {
			search_K_multi(multi_K);
			if(eval_best > multi_K[obj*10]) {
				accept = 1;
				eval_best = multi_K[obj*10];
			}
		}
		if(!accept) {
			// 元に戻す．
			#ifdef REPLACE_WHOLE
				memcpy(pivot.piv0[dim], piv_temp, db_header->element_size * db_header->data_dim);
				pivot.rad[dim] = rad_temp;
			#else
				for(int j = 0; j < f_num; j++) {
					pivot.piv0[dim][axi[j]] = 255 - pivot.piv0[dim][axi[j]];	// 0 <=> 255 のフリップ
				}
				pivot.rad[dim] = compute_rad(pivot.piv0[dim], db_header->data_dim, MAX_SAMPLESIZE / 10, count_tbl);  // 半径計算 (BP と同じ計算方法の場合)
			#endif
			temp_sketch = sketch; sketch = y_sketch; y_sketch = temp_sketch; // 元に戻す．
			remake_query_sketch_for_resample(dim, resampling_size);
		} else {
			num_accept++;
		}
		if(final_check && not_improved > PJT_DIM * FTR_DIM) {
			fprintf(stderr, "Convergence at i = %d\n", i - PJT_DIM * FTR_DIM);
			printf("Convergence at i = %d\n", i - PJT_DIM * FTR_DIM);
			break;
		}
	}
	
	fflush(stdout);
	free(ans);
	free(sketch);
	free(y_sketch);
	sketch = y_sketch = NULL;
	fprintf(stderr, "\nNum_accept = %6d\n", num_accept);
	fprintf(stderr, "\nFinal evaluation = %6.2lf\n", eval_best);
	printf("\nFinal evaluation = %6.2lf\n", eval_best);
}

unsigned int data_to_sketch(ftr_type o, int dim)
{
	int j;
	sketch_type sk = 0;
	for(j = 0; j < dim; j++) {
		#if SKETCH_PARTITION != 3
			write_bit(j, DISTANCE_2(o, pivot.piv0[j], db_header->data_dim, pivot.rad[j]) <= pivot.rad[j], &sk); 
		#else
			write_bit(j, pdist_L1(o, pivot.piv0[j], j * 4, (j + 1) * 4 - 1) <= pivot.rad[j], &sk); 
		#endif
	}
	return sk;
}

void data_to_sketch_1bit(ftr_type o, int dim, int idx)
{
	#if SKETCH_PARTITION != 3
		write_bit(dim, DISTANCE(o, pivot.piv0[dim], db_header->data_dim)  <= pivot.rad[dim], &sample_sketch[idx]);
	#else
		write_bit(dim, pdist_L1(o, pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &sample_sketch[idx]);
	#endif
}

void make_sketch(void)
{
	int i, j;
	if(sketch == NULL) 
		sketch = (sketch_type *)malloc(sizeof(sketch_type) * db_header->data_num);
	for(i = 0; i < db_header->data_num; i++) {
		sketch[i] = 0;
	}

	for(j = 0; j < PJT_DIM; j++) {
		set_dist_L2_22(pivot.piv0[j]);
		#pragma omp parallel for
		for(i = 0; i < db_header->data_num; i++) {
			#if SKETCH_PARTITION != 3
				write_bit(j, dist_L2_22(database[i])  <= pivot.rad[j], &sketch[i]);
			#else
				write_bit(j, pdist_L1(database[i], pivot.piv0[j], j * 4, (j + 1) * 4 - 1)  <= pivot.rad[j], &sketch[i]);
			#endif
		}
	}
}


// sketch[i] を temp_sketch[i] にコピーして，dim で指定したビットだけ書き直す．
void make_y_sketch(int dim)
{
	int i;

	set_dist_L2_22(pivot.piv0[dim]);
	#pragma omp parallel for
	for(i = 0; i < db_header->data_num; i++) {
		y_sketch[i] = sketch[i];
		#if SKETCH_PARTITION != 3
			write_bit(dim, dist_L2_22(database[i])  <= pivot.rad[dim], &y_sketch[i]);
		#else
			write_bit(dim, pdist_L1(database[i], pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &y_sketch[i]);
		#endif
	}
}

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

void make_query_sketch(void)
{
	int i;
	#pragma omp parallel for
	for(i = 0; i < query_num; i++) {
		query[i].sketch = data_to_sketch(query[i].data, PJT_DIM);
	}
}

void make_query_sketch_for_resample(int resampling_size)
{
	int i;
	int j;
	int l;

	#ifdef EVAL_BY_SEQUENTIAL_FILTERING
	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++) {
		for(int p = 0; p < 4; p++) {
			for(int n = 0; n < 256; n++) {
				query[i].tbl[p][n] = 0;
			}
		}
		query[i].answer_sketch = 0;
	}
	#endif

	for(j = 0; j < PJT_DIM; j++) {
		set_dist_L2_22(pivot.piv0[j]);
		#pragma omp parallel for
		for(i = 0; i < resampling_size; i++) {
			#if SKETCH_PARTITION != 3
				dist_type dist = dist_L2_22(query[i].data);
				write_bit(j, dist <= pivot.rad[j], &(query[i].sketch));
			#else
				fprintf(stderr, "Not implemented for SKETCH_PARTITON = 3\n");
				exit(0);
				//write_bit(j, pdist_L1(database[i], pivot.piv0[j], j * 4, (j + 1) * 4 - 1)  <= pivot.rad[j], &sketch[i]);
			#endif
			query[i].bd[j] = abs(dist - pivot.rad[j]);
			for(l = j - 1; l >= 0 && query[i].bd[query[i].idx[l]] > query[i].bd[j]; l--) {
				query[i].idx[l + 1] = query[i].idx[l];
			}
			query[i].idx[l + 1] = j;
			#ifdef EVAL_BY_SEQUENTIAL_FILTERING
			int tn = j / 8, bp = j % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
			for(int n = 0; n < 256; n++) {
				if(n & (1 << bp)) {
					query[i].tbl[tn][n] += query[i].bd[j];
				}
			}
			#if SKETCH_PARTITION != 3
				dist_type ans_dist = dist_L2_22(query[i].answer_data);
				write_bit(j, ans_dist <= pivot.rad[j], &(query[i].answer_sketch));
			#else
				fprintf(stderr, "Not implemented for SKETCH_PARTITON = 3\n");
				exit(0);
				//write_bit(j, pdist_L1(database[i], pivot.piv0[j], j * 4, (j + 1) * 4 - 1)  <= pivot.rad[j], &sketch[i]);
			#endif
			#endif
		}
	}
	#ifdef EVAL_BY_SEQUENTIAL_FILTERING
	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++) {
		query[i].answer_score = priority(query[i].answer_sketch, &query[i]);
	}
	#endif
//	if(!first_flag && resampling_size == sample_size) {
//		printf("Prepared All Sketches!\n");
//		prepared_all_sketches = 1;
//	}
//	first_flag = 0;
}

void remake_query_sketch(int dim)
{
	int i;

	#pragma omp parallel for
	for(i = 0; i < query_num; i++) {
#if SKETCH_PARTITION != 3
		write_bit(dim, DISTANCE(query[i].data, pivot.piv0[dim], db_header->data_dim)  <= pivot.rad[dim], &query[i].sketch);
#else
		write_bit(dim, pdist_L1(query[i].data, pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &query[i].sketch);
#endif
	}
}

void remake_query_sketch_for_resample(int dim, int resampling_size)
{
	int i;

//	if(prepared_all_sketches) {
//		resampling_size = sample_size;
//	}
	set_dist_L2_22(pivot.piv0[dim]);
	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++) {
		dist_type dist;
		int l;

		#if SKETCH_PARTITION != 3
			dist = dist_L2_22(query[i].data);
			write_bit(dim, dist <= pivot.rad[dim], &(query[i].sketch));
		#else
			fprintf(stderr, "Not implemented for SKETCH_PARTITON = 3\n");
			exit(0);
			// write_bit(dim, pdist_L1(query[i].data, pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &query[i].sketch);
		#endif
		#ifdef EVAL_BY_SEQUENTIAL_FILTERING
		dist_type bd_old = query[i].bd[dim];
		dist_type bd_new = abs(dist - pivot.rad[dim]);
		int tn = dim / 8, bp = dim % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
		for(int n = 0; n < 256; n++) {
			if(n & (1 << bp)) {
				query[i].tbl[tn][n] += (bd_new - bd_old);
			}
		}
		query[i].bd[dim] = bd_new;
		#if SKETCH_PARTITION != 3
			dist = dist_L2_22(query[i].answer_data);
			write_bit(dim, dist <= pivot.rad[dim], &(query[i].answer_sketch));
			query[i].answer_score = priority(query[i].answer_sketch, &query[i]);
		#else
			fprintf(stderr, "Not implemented for SKETCH_PARTITON = 3\n");
			exit(0);
			// write_bit(dim, pdist_L1(query[i].data, pivot.piv0[dim], dim * 4, (dim + 1) * 4 - 1)  <= pivot.rad[dim], &query[i].sketch);
		#endif
		#else
		query[i].bd[dim] = abs(dist - pivot.rad[dim]);
		#endif

		for(l = 0; l < PJT_DIM; l++) {
			if(query[i].idx[l] == dim) {
				if(l == 0 || query[i].bd[query[i].idx[l - 1]] < query[i].bd[dim]) { // 右に挿入
					for( ; l < PJT_DIM - 1 && query[i].bd[query[i].idx[l + 1]] < query[i].bd[dim]; l++) {
						query[i].idx[l] = query[i].idx[l + 1];
					}
					query[i].idx[l] = dim;
				} else { // 左に挿入
					for( ; l > 0 && query[i].bd[query[i].idx[l - 1]] > query[i].bd[dim]; l--) {
						query[i].idx[l] = query[i].idx[l - 1];
					}
					query[i].idx[l] = dim;
				}
				break;
			}
		}
	}
}

int bit_count(sketch_type a)
{
	int c = 0;

	if(sizeof(sketch_type) == sizeof(short))
		c = bitcnt_tbl[a]; // 16bit
	else {
		//for(; a; a >>= 8)  // 64bit等にも対応　ただし遅い
		//	c += bitcnt_tbl[a & 0xff]; 
		c = bitcnt_tbl[a & 0xff] + bitcnt_tbl[(a >> 8) & 0xff] + bitcnt_tbl[(a >> 16) & 0xff] + bitcnt_tbl[(a >> 24) & 0xff]; // 32bit
	}
	
	return c;
}

#if SKETCH_PARTITION != 3

void bound(ftr_type q, int b[])
{
	int j;
	
	for(j = 0; j < PJT_DIM; j++) {
		//b[j] = abs(DISTANCE(pivot.piv0[j], q, db_header->data_dim) - DISTANCE(pivot.piv0[j], med, db_header->data_dim));
		b[j] = abs(DISTANCE(pivot.piv0[j], q, db_header->data_dim) - pivot.rad[j]);
#ifdef EUCLID
//		fprintf(stderr, "EUCLID"); exit(0);
		if(lb == 10) {
			b[j] = sqrt(b[j]);
		}
#endif
	}
}

void bound_with_idx(ftr_type q, bd_type b[])
{
	int j;
	
	for(j = 0; j < PJT_DIM; j++) {
		b[j].idx = j;
		//b[j] = abs(DISTANCE(pivot.piv0[j], q, db_header->data_dim) - DISTANCE(pivot.piv0[j], med, db_header->data_dim));
		b[j].bd = abs(DISTANCE(pivot.piv0[j], q, db_header->data_dim) - pivot.rad[j]);
#ifdef EUCLID
//		fprintf(stderr, "EUCLID"); exit(0);
		if(lb == 10) {
			b[j].bd = sqrt(b[j].bd);
		}
#endif
	}
}
#else
void bound(ftr_type q, int b[])
{
	int j;
	
	for(j = 0; j < PJT_DIM; j++) {
		//b[j] = abs(DISTANCE(pivot.piv0[j], q, db_header->data_dim) - DISTANCE(pivot.piv0[j], med, db_header->data_dim));
		b[j] = abs(pdist_L1(pivot.piv0[j], q, j * 4, (j + 1) * 4 - 1) - pivot.rad[j]);
	}
}
#endif

double collision(double *collsize, int idx[], int samplesize) {
	int i;
	double n;       // 衝突している組の個数
 	double score = 0;       // スコア
	int p = 0;
	int temp[MAX_SAMPLESIZE];

	for(i = 0; i < samplesize; i++)
		temp[i] = ind_re[i];
	qsort(temp, samplesize, sizeof(int), comp_sketch); // sampleのindexをsketchでソート

	for (i = 0; i < samplesize; i++) {
		n = 0;
		while(i < samplesize - 1 && sample_sketch[temp[i]] == sample_sketch[temp[i + 1]]) {
			idx[p] = temp[i];
			n++;
			p++;
			i++;
		}
		if(n > 0) {
			idx[p] = temp[i];
			n++;
			p++;
			score += n * (n - 1) / 2;
		}
	}
	*collsize = p;
	return score;
}

double collision2(double coll_size, int coll_idx[], int piv) {
	int i;
#if SKETCH_PARTITION == 3 
	int dim_from = piv * 4, dim_to = (piv + 1) * 4 - 1;
#endif
	double n0, n1;
	double sc;
	sc = 0;
	for(i = 0; i < coll_size; i++){ 
		n0 = 0;
		n1 = 0;
		while(i < coll_size - 1 && sample_sketch[coll_idx[i]] == sample_sketch[coll_idx[i + 1]]){
			if(flag == 1){
#if SKETCH_PARTITION != 3 
				if(DISTANCE(sample[coll_idx[i]], pivot.piv0[piv], db_header->data_dim) <= pivot.rad[piv]) {
#else
				if(pdist_L1(sample[coll_idx[i]], pivot.piv0[piv], dim_from, dim_to) <= pivot.rad[piv]) {
#endif
					n0++;
				}
				else {
					n1++;
				}
			}
			if(flag == 0){
				if(b_dist[coll_idx[i]][piv] <= pivot.rad[piv]) {
					n0++;
				}
				else {
					n1++;
				}
			}
		   
		i++;
		}

		if(flag == 1){
#if SKETCH_PARTITION != 3 
			if(DISTANCE(sample[coll_idx[i]], pivot.piv0[piv], db_header->data_dim) <= pivot.rad[piv]) {
#else
			if(pdist_L1(sample[coll_idx[i]], pivot.piv0[piv], dim_from, dim_to) <= pivot.rad[piv]) {
#endif
				n0++;
			}
			else {
				n1++;
			}
		}
		
		if(flag == 0){
			if(b_dist[coll_idx[i]][piv] <= pivot.rad[piv]) {
				n0++;
			}
			else {
				n1++;
			}
		}
		
		
		sc += (n0*(n0-1) + n1*(n1-1)) / 2;
	}

	return sc;
}

double correl(unsigned int a[], unsigned int b[], int n) 
{
	int i;
	double sa, sb, sa2, sb2, ava, avb, va, vb, cov;
//	sum ( x - a )^2 = sum (x^2) - 2* sum(x * a) + n * a^2
//                  = sum x^2 - 2n a^2 + na^2 = sum x^2 - na^2
	
	sa = sb = sa2 = sb2 = 0;
	for(i = 0; i < n; i++) {
		sa += a[i];
		sa2 += a[i] * a[i];
		sb += b[i];
		sb2 += b[i] * b[i];
	}
	ava = sa / n;
	va = sa2 - n * ava * ava;
	avb = sb / n;
	vb = sb2 - n * avb * avb;
	cov = 0;
	for(i = 0; i < n; i++) {
		cov += (a[i] - ava) * (b[i] - avb);
	}

//	printf("ava = %8.6f, avb = %8.6f, sa2 = %10.2f, sb2 = %10.2f", ava, avb, sa2, sb2);
	return cov / sqrt(va * vb);
	
}

int comp_sketch(const void *a, const void *b) {

	if(sample_sketch[*(int *)a] < sample_sketch[*(int *)b])
		return -1;
	else if(sample_sketch[*(int *)a] == sample_sketch[*(int *)b])
		return 0;
	else
		return 1;
	
}
/*
int comp_query(const void *a, const void *b) {

	if(a < b)
		return 1;
	else if(a == b)
		return 0;
	else
		return -1;
	
}
*/
int cmp_buff( const void *p, const void *q ) {
	if(((buff_type*)p)->dist < ((buff_type*)q)->dist)
		return -1;
	else if(((buff_type*)p)->dist < ((buff_type*)q)->dist)
		return 0;
	else return 1;

}

/*
void make_table2(int b[], int table[][BIT]) {
	int i, j, k;
	int bit;

	for(i = 0; i < BIT; i++) { // 0 ~ 255のテーブル作成
		for(j = 0; j < 8; j++) { // 8bit
			for(k = 0; k < 4; k++) {
				bit = j + k * 8;
				if((i >> j) & 1) {
					if(lb == 3) {
						if(table[k][i] < b[bit])
							table[k][i] = b[bit];
					}
					else if(lb == 1) {
						table[k][i] += b[bit];
					}
					else {
						table[k][i] += b[bit] * b[bit];
					}
				}
			}
		}
	}
}

void make_table3(int b_in[], int b_out[], int table[][BIT], sketch_type q_sk) {
	int i, j, k;
	int bit;

	for(i = 0; i < BIT; i++) { // 0 ~ 255のテーブル作成
		for(j = 0; j < 8; j++) { // 8bit
			for(k = 0; k < 4; k++) {
				bit = j + k * 8;
				if((i >> j) & 1) { // sketchが衝突しない場合
					if((q_sk >> bit) & 1) { // queryが内側
						if(lb == 3) { // Linf
							if(table[k][i] < b_out[bit])
								table[k][i] = b_out[bit];
						}

						else if(lb == 1) { // L1
							table[k][i] += b_out[bit];
						}
					}
					else {				  // queryが外側
						if(lb == 3) { // Linf
							if(table[k][i] < b_in[bit])
								table[k][i] = b_in[bit];
						}
						else if(lb == 1) { // L1
							table[k][i] += b_in[bit];
						}
					}
				}
				
				else { // sketchが衝突する場合
					if((q_sk >> bit) & 1) { // queryが内側
						if(lb == 3) { // Linf
							if(table[k][i] < b_in[bit])
								table[k][i] = b_in[bit];
						}
						else if(lb == 1) { // L1
							table[k][i] += b_in[bit];
						}
					}
					else {				  // queryが外側
						if(lb == 3) { // Linf
							if(table[k][i] < b_out[bit])
								table[k][i] = b_out[bit];
						}
						else if(lb == 1) { // L1
							table[k][i] += b_out[bit];
						}
					}
				}
				
			}
		}
	}
}


dist_type search_score(sketch_type s, table_type table) {

	if (lb == 0) // Hamming
		return bit_count(s);
	else {
		if (lb == 3) { // L_inf
			int temp0, temp1;
			temp0 = (table.t0[s & 0xff] > table.t1[(s >> 8) & 0xff]) ? table.t0[s & 0xff] : table.t1[(s >> 8) & 0xff];
			temp1 = (table.t2[(s >> 16) & 0xff] > table.t3[(s >> 24) & 0xff]) ? table.t2[(s >> 16) & 0xff] : table.t3[(s >> 24) & 0xff];
			return (temp0 > temp1) ? temp0 : temp1;
		}
		else // L1, L2
			return table.t0[s & 0xff] + table.t1[(s >> 8) & 0xff] + table.t2[(s >> 16) & 0xff] + table.t3[(s >> 24) & 0xff];
	}

}
*/

void make_ham_table(void) {
	int i;
	unsigned int v;

	for(i = 0; i < BIT; i++) {
		v = i;
		for (bitcnt_tbl[i] = 0; v; v >>= 1) {
			bitcnt_tbl[i] += v & 1;
		}
	}
}

void ham_search_16bit(int search_mode) {
	// search_mode = 0  通常の検索（ソートしたデータベースを用いる）
	// search_mode = 1  通常の検索（ソートしないデータベースを用いる）
	// search_mode = 2  通常の検索（スケッチの列挙のみで実際の検索をしない）
		
	int i, j, q, k;
	int *idx;
	idx = (int *)malloc(sizeof(int) * db_header->data_num);
	int *sk_st, *sk_num;
	sk_st = (int *)calloc(BIT ,sizeof(int));
	sk_num = (int *)calloc(BIT, sizeof(int));
	dist_type dist = 0;
	int *bit_idx;
	bit_idx = (int *)malloc(BIT * sizeof(int));
	sketch_type s;
#if SKETCH_PARTITION == 3
// 安全に枝刈りされた回数を記録
	int pruned = 0, total_pruned = 0, skipped = 0, total_skipped = 0, hamm = 0;
#endif

	for(i = 0; i < db_header->data_num; i++) {
		idx[i] = i;
		sk_num[sketch[i]]++;
	}

	qsort(idx, db_header->data_num, sizeof(int), comp_sk); // sketchでindexをソート
	
	if(database2 == NULL) {
		alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
	}
	for(i = 0; i < db_header->data_num; i++) {
		memcpy(database2[i], database[idx[i]], FTR_DIM);
	}
	
	int c = 0;
	for(i = 0; i < BIT; i++) {
		sk_st[i] = c;  // ソートされたsketchのそれぞれの開始位置を記憶
		c += sk_num[i];
		bit_idx[i] = i ;
	}

	qsort(bit_idx, BIT, sizeof(int), comp_bit);
	
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;
	clock_gettime(CLOCK_REALTIME, &tp1);
	
	sk_c = sk_c_0 = 0;
	for(q = 0; q < query_num; q++) {
		fprintf(stderr, "query = %5d\r", q + 1);	
#if SKETCH_PARTITION == 3
		pruned = skipped = hamm = 0;
#endif
		for(i = 0, k = 0; i < BIT && k < K; i++) {
			s = query[q].sketch ^ bit_idx[i];
#if SKETCH_PARTITION == 3
			if(hamm < bit_count(bit_idx[i])) {
				hamm = bit_count(bit_idx[i]);
				//fprintf(stderr, "hamm = %d, hamm * expand = %d, dist = %d\n", hamm, hamm * (hamm <= P_DIST_FAR ? EXPAND_1 : EXPAND_2), ans[q].dist);
				if(SKIP_CONDITION(hamm, ans[q].dist)) {
					skipped++;
					// fprintf(stderr, "skipped for q = %d, at k = %d, Hamm = %d, ans.dist = %d\n",q + 1, k, hamm, ans[q].dist);
					break;
				}
			}
#endif
			sk_c_0 += (sk_num[s] == 0);
			for(j = 0; j < sk_num[s] && k < K; j++, k++) {
				if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[sk_st[s] + j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[sk_st[s] + j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				if(dist < ans[q].dist) {
					ans[q].dist = dist;
					ans[q].idx = idx[sk_st[s] + j];
				}
			}
		}
		sk_c += i;
#if SKETCH_PARTITION == 3
		total_pruned += pruned;
		total_skipped += skipped;
#endif
	}
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
#if SKETCH_PARTITION == 3
	fprintf(stderr,"total_pruned = %d, total_skipped = %d\n", total_pruned, total_skipped);
#endif
}

void Linf_search_16bit(int search_mode) {
	int i, j, q, k;
	int *idx;
	idx = (int *)malloc(sizeof(int) * db_header->data_num);
	int *sk_st, *sk_num;
	sk_st = (int *)calloc(BIT ,sizeof(int));
	sk_num = (int *)calloc(BIT, sizeof(int));
	int b_idx[PJT_DIM];
	dist_type dist = 0;
	sketch_type s;

	for(i = 0; i < db_header->data_num; i++) {
		idx[i] = i;
		sk_num[sketch[i]]++;
	}

	qsort(idx, db_header->data_num, sizeof(int), comp_sk); // sketchでindexをソート
	
	int c = 0;
	for(i = 0; i < BIT; i++) {
		sk_st[i] = c;  // ソートされたsketchのそれぞれの開始位置を記憶
		c += sk_num[i];
	}
	
	if(database2 == NULL) {
		alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
	}
	for(i = 0; i < db_header->data_num; i++) {
		memcpy(database2[i], database[idx[i]], FTR_DIM);
	}

	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;

	sk_c = sk_c_0 =0;
	clock_gettime(CLOCK_REALTIME, &tp1);
	for(q = 0; q < query_num; q++) {
		// fprintf(stderr, "query = %5d\r", q + 1);
		bound(query[q].data, bd);
		
		for(i = 0; i < PJT_DIM; i++)
			b_idx[i] = i;
		qsort(b_idx, PJT_DIM, sizeof(int), comp_b);
		
		s = query[q].sketch;
		k = 0;
		
		for(j = 0; j < sk_num[s] && k < K; j++, k++) {
			if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[sk_st[s] + j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[sk_st[s] + j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[sk_st[s] + j];
			}
			
		}
		
		for(i = 0; i < BIT && k < K; i++) {
			s ^= 1 << b_idx[bitcnt_tbl[i ^ (i + 1)] - 1];	
			sk_c_0 += (sk_num[s] == 0);
			for(j = 0; j < sk_num[s] && k < K; j++, k++) {
				if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[sk_st[s] + j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[sk_st[s] + j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				if(dist < ans[q].dist) {
					ans[q].dist = dist;
					ans[q].idx = idx[sk_st[s] + j];
				}
				
			}
		}
		
		sk_c += i;
		
	}
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
}

// lb = 1 				score_1
// lb = 10,11, 12, ...	score_p (p = lb/10.0)
void L1_search_16bit_2(int search_mode) { // 列挙による高速化あり
	int i, q;
	int skipped = 0, total_skipped = 0;
	
	prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);
	alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
	sort_database2();

	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;
	sk_c = sk_c_0 =0;
	clock_gettime(CLOCK_REALTIME, &tp1);

	for(q = 0; q < query_num; q++) {
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE qu, qu2;
		struct_que *que;
		dist_type dist = 0;
		int i, k;

		que = (struct_que *)malloc(sizeof(struct_que));
		bound_with_idx(query[q].data, bd_idx);
		if(lb >= 5) {
			p = (double)lb / 10;
			for(i = 0; i < PJT_DIM; i++)
				bd_idx[i].bd = pow(bd_idx[i].bd, p);
		}
		qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);
		s = query[q].sketch;
// printf("s = %d = ", s); print_bin(s); printf("\n");
		k = 0;
		skipped = 0;
		for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		}
		i = 1;
		que->qsize = 0;
		qu.key = bd_idx[0].bd;
		qu.sk = s ^ (1 << bd_idx[0].idx);
		qu.idx = 1;
		ENQ_P(&qu, que);

		while(DEQ_P(&qu, que) && k < K) {
// printf("s = %d = ", qu.sk); print_bin(qu.sk); printf("\n");
			sk_c_0 += (bkt[qu.sk + 1] - bkt[qu.sk] == 0);
			for(int j = bkt[qu.sk]; j < bkt[qu.sk + 1] && k < K; j++, k++) {
				if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				if(dist < ans[q].dist) {
					ans[q].dist = dist;
					ans[q].idx = idx[j];
				}
			}
			if(qu.idx < PJT_DIM) {
				s = qu.sk ^ (1 << bd_idx[qu.idx].idx);
				qu2.key = qu.key + bd_idx[qu.idx].bd;
				qu2.sk = s;
				qu2.idx = qu.idx + 1;
#ifdef SKIP_ENQ
				if(SKIP_CONDITION(qu2.key, ans[q].dist)) {
					skipped++;
				} else {
					ENQ_P(&qu2, que);
				}
#else
				ENQ_P(&qu2, que);
#endif
				qu2.key = qu.key + bd_idx[qu.idx].bd - bd_idx[qu.idx - 1].bd;
				qu2.sk = s ^ (1 << bd_idx[qu.idx - 1].idx);
				qu2.idx = qu.idx + 1;
#ifdef SKIP_ENQ
				if(SKIP_CONDITION(qu2.key, ans[q].dist)) {
					skipped++;
				} else {
					ENQ_P(&qu2, que);
				}
#else
				ENQ_P(&qu2, que);
#endif
			}
			i++;
		}
		sk_c += i;
		free(que);
		total_skipped += skipped;
	}
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
	fprintf(stderr, "total_skipped = %d\n", total_skipped);
//	free_database2(db_header->data_num);
}

/*
void debug(int idx, ftr_type q, ftr_type d) {
	printf("idx = %10d\n", idx);
	for(int i = 0; i < db_header->data_dim; i += 20) {
		printf("(q) dim = %4d: ", i);
		for(int j = i; j < i + 20; j++) {
			printf("%4d", q[j]);
		}
		printf("\n");
	}
	for(int i = 0; i < db_header->data_dim; i += 20) {
		printf("(d) dim = %4d: ", i);
		for(int j = i; j < i + 20; j++) {
			printf("%4d", d[j]);
		}
		printf("\n");
	}
	printf("dist = %d\n", DISTANCE_2(q, d, db_header->data_dim, 100000000));
	fflush(stdout);
	exit(0);
}
*/

void L1_search_16bit_2_c2_n(int search_mode) { // 改善版の列挙による高速化あり
	// search_mode = 0  通常の検索（ソートしたデータベースを用いる）
	// search_mode = 1  通常の検索（ソートしないデータベースを用いる）
	// search_mode = 2  スケッチの列挙のみ（実際の検索をしない）
	// search_mode = 3  通常の検索（配列を一つでソートしたデータベースを用いる）ここでのみ
	// search_mode = 4  通常の検索（HDDに保存したbucketデータベースを用いる）ここでのみ (オンメモリ）
	// search_mode = 5  通常の検索（HDDに保存したbucketデータベースを用いる）ここでのみ（ftr は2次記憶上のまま）
	int i, q;
	ftr_type ftr_data =  (ftr_type)malloc(db_header->element_size * db_header->data_dim);
	sub_fb_type data_id;
	struct timespec tp3, tp4;

	if(search_mode == 5) {
		// 特徴データ全体をdatabaseに読み込まないで，照合のたびに読み込む
//		ftr_data = (ftr_type)malloc(db_header->element_size * db_header->data_dim);
	} else if(search_mode == 4) {
		search_mode = 3;
	} else {
		prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);
		if(search_mode == 0) {
			alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
			sort_database2();
		} else if(search_mode == 3) {
			sort_database_self();
		}
	}
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;
	sk_c = sk_c_0 =0;
	clock_gettime(CLOCK_REALTIME, &tp1);

	for(q = 0; q < query_num; q++) {
//	for(q = 0; q < 1; q++) {
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE_c2 qu, qu2;
		struct_que_c2_n *que = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n));
		dist_type dist = 0;
		int i = 0, k;

		clock_gettime(CLOCK_REALTIME, &tp3);
		bound_with_idx(query[q].data, bd_idx);
		if(lb >= 5) {
			p = (double)lb / 10;
			for(i = 0; i < PJT_DIM; i++)
				bd_idx[i].bd = pow(bd_idx[i].bd, p);
		}
		qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

		k = 0;

		s = query[q].sketch;
// printf("s = %d = ", s); print_bin(s); printf("\n");
		for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			if(search_mode == 5) {
				if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
					fprintf(stderr, "get_ftr error at %d\n", j);
					exit(0);
				}
				if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
					fprintf(stderr, "cannot read sub_fb at %d\n", j);
					exit(0);
				}
				dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 3) {
				dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		}
		
		s = s ^ (1 <<  bd_idx[0].idx);
// printf("s = %d = ", s); print_bin(s); printf("\n");
		for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			if(search_mode == 5) {
				if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
					fprintf(stderr, "get_ftr error at %d\n", j);
					exit(0);
				}
				if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
					fprintf(stderr, "cannot read sub_fb at %d\n", j);
					exit(0);
				}
				dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 3) {
//				if(q == 0 && idx[sk_st[s] + j] == 683296) debug(idx[sk_st[s] + j], database[sk_st[s] + j], query[q].data);
				dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		}

		make_empty_que_c2_n(que);

		// enq pattern of 0...10
		qu.cursor = new_que_e2_n(que);
		qu.key = bd_idx[1].bd;
		que->details[qu.cursor].sk = query[q].sketch ^ (1 << bd_idx[1].idx);
		que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
		enq_c2_n(&qu, que);		

		while(deq_c2_n(&qu, que) && k < K) {
			s = que->details[qu.cursor].sk;
// printf("s = %d = ", s); print_bin(s); printf("\n");
			sk_c_0 += (bkt[s] == bkt[s + 1]);
			for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
				if(search_mode == 5) {
					if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
						fprintf(stderr, "get_ftr error at %d\n", j);
						exit(0);
					}
					if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
						fprintf(stderr, "cannot read sub_fb at %d\n", j);
						exit(0);
					}
					dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
				} else if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else if(search_mode == 3) {
//					if(q == 0 && idx[sk_st[s] + j] == 683296) debug(idx[sk_st[s] + j], database[sk_st[s] + j], query[q].data);
					dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				if(dist < ans[q].dist) {
					ans[q].dist = dist;
					ans[q].idx = idx[j];
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
						qu2.key = qu.key + bd_idx[m + 1].bd - bd_idx[m].bd;
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[m + 1].idx)) ^ (1 << bd_idx[m].idx);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
						// Y010^m -> Y010^{m-1}1
						qu.key = qu.key + bd_idx[0].bd;
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
						enq_c2_n(&qu2, que);
					} else {
						qu.key = qu.key + bd_idx[0].bd;
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
					}
				}
				break;
			case 4:  // X0100 -> enq(X0101) and enq(X1000)
				// X1000
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + bd_idx[3].bd - bd_idx[2].bd;
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[3].idx)) ^ (1 << bd_idx[2].idx);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
				// X0101
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 1:  // X0001 -> enq(X0010)
			case 5:  // X0101 -> enq(X0110)
			case 9:  // X1001 -> enq(X1010)
			case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
				qu.key = qu.key + bd_idx[1].bd - bd_idx[0].bd;
				que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[1].idx)) ^ (1 << bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 2:  // X0010 -> enq(X0011) and enq(X0100)
			case 10: // X1010 -> enq(X1011) and enq(X1100)
				// X0100 and X1100
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key +  bd_idx[2].bd -  bd_idx[1].bd;
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[2].idx)) ^ (1 << bd_idx[1].idx);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
				// X0011 and X1011
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 6:  // X0110 -> enq(X0111)
			case 12: // X1100 -> enq(X1101)
			case 14: // X1110 -> enq(10111)
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 3:  // X0011
			case 7:  // X0111
			case 11: // X1011
			case 15: // X1111 -> nothing to do
				break;
			}
			i++;
		}
		sk_c += i;
		free(que);
		clock_gettime(CLOCK_REALTIME,&tp4);
		sec = tp4.tv_sec - tp3.tv_sec;
		nsec = tp4.tv_nsec - tp3.tv_nsec;
		if(nsec < 0){
			sec--;
			nsec += 1000000000L;
		}
	#ifndef NO_INTERMEDIATE_RESULT
		printf("q = %d, dist = %d, idx = %d, %ld.%09ld (sec)\n", q, ans[q].dist, ans[q].idx, sec, nsec);
	#endif
	}
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
}

void L1_search_16bit_2_c2_n_new_bucket(int search_mode) { // 改善版の列挙による高速化あり
	// search_mode = 0  通常の検索（ソートしたデータベースを用いる）
	// search_mode = 1  通常の検索（ソートしないデータベースを用いる）
	// search_mode = 2  スケッチの列挙のみ（実際の検索をしない）
	// search_mode = 3  通常の検索（配列を一つでソートしたデータベースを用いる）ここでのみ
	// search_mode = 4  通常の検索（HDDに保存したbucketデータベースを用いる）ここでのみ (オンメモリ）
	// search_mode = 5  通常の検索（HDDに保存したbucketデータベースを用いる）ここでのみ（ftr は2次記憶上のまま）
	int i, q;
	ftr_type ftr_data =  (ftr_type)malloc(db_header->element_size * db_header->data_dim);
	sub_fb_type data_id;
	struct timespec tp3, tp4;

	if(search_mode == 5) {
		// 特徴データ全体をdatabaseに読み込まないで，照合のたびに読み込む
//		ftr_data = (ftr_type)malloc(db_header->element_size * db_header->data_dim);
	} else if(search_mode == 4) {
		search_mode = 3;
	} else {
		prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);
		if(search_mode == 0) {
			alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
			sort_database2();
		} else if(search_mode == 3) {
			sort_database_self();
		}
	}
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;
	sk_c = sk_c_0 =0;
	clock_gettime(CLOCK_REALTIME, &tp1);

#ifdef Q_NUM_FROM
	int query_from = Q_NUM_FROM;
#else
	int query_from = 0;
#endif
	for(q = query_from; q < query_from + query_num; q++) {
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE_c2 qu, qu2;
		struct_que_c2_n *que = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n));
		dist_type dist = 0;
		int i = 0, j, k;

		clock_gettime(CLOCK_REALTIME, &tp3);
		bound_with_idx(query[q].data, bd_idx);
		if(lb >= 5) {
			p = (double)lb / 10;
			for(i = 0; i < PJT_DIM; i++)
				bd_idx[i].bd = pow(bd_idx[i].bd, p);
		}
		qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

		k = 0;

		s = query[q].sketch;
		for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			if(search_mode == 5) {
				if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
					fprintf(stderr, "get_ftr error at %d\n", j);
					exit(0);
				}
				if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
					fprintf(stderr, "cannot read sub_fb at %d\n", j);
					exit(0);
				}
				dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 3) {
				dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		}
		
		s = s ^ (1 <<  bd_idx[0].idx);
		for(j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			if(search_mode == 5) {
				if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
					fprintf(stderr, "get_ftr error at %d\n", j);
					exit(0);
				}
				if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
					fprintf(stderr, "cannot read sub_fb at %d\n", j);
					exit(0);
				}
				dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 3) {
				dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		}

		make_empty_que_c2_n(que);

		// enq pattern of 0...10
		qu.cursor = new_que_e2_n(que);
		qu.key = bd_idx[1].bd;
		que->details[qu.cursor].sk = query[q].sketch ^ (1 << bd_idx[1].idx);
		que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
		enq_c2_n(&qu, que);		

		while(deq_c2_n(&qu, que) && k < K) {
			s = que->details[qu.cursor].sk;
			sk_c_0 += (bkt[s] == bkt[s + 1]);
			for(j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
				if(search_mode == 5) {
					if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
						fprintf(stderr, "get_ftr error at %d\n", j);
						exit(0);
					}
					if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
						fprintf(stderr, "cannot read sub_fb at %d\n", j);
						exit(0);
					}
					dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
				} else if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else if(search_mode == 3) {
					dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				if(dist < ans[q].dist) {
					ans[q].dist = dist;
					ans[q].idx = idx[j];
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
						qu2.key = qu.key + bd_idx[m + 1].bd - bd_idx[m].bd;
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[m + 1].idx)) ^ (1 << bd_idx[m].idx);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
						// Y010^m -> Y010^{m-1}1
						qu.key = qu.key + bd_idx[0].bd;
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
						enq_c2_n(&qu2, que);
					} else {
						qu.key = qu.key + bd_idx[0].bd;
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
					}
				}
				break;
			case 4:  // X0100 -> enq(X0101) and enq(X1000)
				// X1000
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + bd_idx[3].bd - bd_idx[2].bd;
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[3].idx)) ^ (1 << bd_idx[2].idx);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
				// X0101
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 1:  // X0001 -> enq(X0010)
			case 5:  // X0101 -> enq(X0110)
			case 9:  // X1001 -> enq(X1010)
			case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
				qu.key = qu.key + bd_idx[1].bd - bd_idx[0].bd;
				que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[1].idx)) ^ (1 << bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 2:  // X0010 -> enq(X0011) and enq(X0100)
			case 10: // X1010 -> enq(X1011) and enq(X1100)
				// X0100 and X1100
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key +  bd_idx[2].bd -  bd_idx[1].bd;
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[2].idx)) ^ (1 << bd_idx[1].idx);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
				// X0011 and X1011
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 6:  // X0110 -> enq(X0111)
			case 12: // X1100 -> enq(X1101)
			case 14: // X1110 -> enq(10111)
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 3:  // X0011
			case 7:  // X0111
			case 11: // X1011
			case 15: // X1111 -> nothing to do
				break;
			}
			i++;
		}
		sk_c += i;
		free(que);
		clock_gettime(CLOCK_REALTIME,&tp4);
		sec = tp4.tv_sec - tp3.tv_sec;
		nsec = tp4.tv_nsec - tp3.tv_nsec;
		if(nsec < 0){
			sec--;
			nsec += 1000000000L;
		}
	#ifndef NO_INTERMEDIATE_RESULT
		printf("q = %d, dist = %d, idx = %d, %ld.%09ld (sec)\n", q, ans[q].dist, ans[q].idx, sec, nsec);
	#endif
	}
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
}

void L1_search_16bit_2_c2_n_new_bucket_para(int search_mode) { // 改善版の列挙による高速化あり，マルチスレッドによる並列検索処理
	// search_mode = 0  通常の検索（ソートしたデータベースを用いる）
	// search_mode = 1  通常の検索（ソートしないデータベースを用いる）
	// search_mode = 2  スケッチの列挙のみ（実際の検索をしない）
	// search_mode = 3  通常の検索（配列を一つでソートしたデータベースを用いる）ここでのみ
	// search_mode = 4  通常の検索（HDDに保存したbucketデータベースを用いる）ここでのみ（オンメモリ）
	// search_mode = 5  通常の検索（HDDに保存したbucketデータベースを用いる）ここでのみ（ftr は2次記憶上のまま）(block_read, parallel)
	// search_mode = 6  スケッチの列挙のみ（HDDに保存したbucketデータベースを用いる）ここでのみ（ftr は2次記憶上のまま）

	fprintf(stderr, "search_mode = %d\n", search_mode);
	
#if ! defined(BLOCK_READ) && ! defined(NUM_THREADS) && defined(SEARCH_ON_SECONDARY)
	ftr_type ftr_data =  (ftr_type)malloc(db_header->element_size * db_header->data_dim);
	sub_fb_type data_id;
#endif

// ブロック単位でまとめ読みして，呼び出し要求に対しては，バッファから次のものを返す
#if defined(BLOCK_READ) || defined(NUM_THREADS)
	// ブロック単位のまとめ読みやマルチスレッドの並列処理のときは，フィルタリング（第1段階検索）と第2段階検索を分離して実行する．
	// BLOCK_READ == 1 の場合でも，スケッチの列挙によるフィルタリング（第1段階検索）と第2段階検索を分離して実行する．
	int *ftr_idx_lst = (int *)malloc(sizeof(int) * K); // 解候補のidx（特徴データのftr内での位置，つまりスケッチ順での順位）のリスト
	int num_idx = 0; // リスト内のデータ数
	#ifndef NUM_THREADS
		fprintf(stderr, "Single thread with block read = %d\n", BLOCK_READ);
		ftr_handle *handle = new_ftr_handle(db_filename, BLOCK_READ);
	#else
		omp_set_num_threads(NUM_THREADS);
		int nt = NUM_THREADS;
		#ifdef BLOCK_READ
		fprintf(stderr, "Parallel, num threads = %d with block read = %d\n", nt, BLOCK_READ);
		ftr_handle *hd_pool[nt];
		for(int t = 0; t < nt; t++) {
			hd_pool[t] = new_ftr_handle(db_filename, BLOCK_READ);
		}
		#else
		fprintf(stderr, "Parallel, num threads = %d without block read\n", nt);
		answer temp_ans[nt];
		#endif
	#endif
#endif

	struct timespec tp3, tp4;

	if(search_mode == 6 || search_mode == 5) {
//		#ifdef BLOCK_READ
//		fprintf(stderr, "BLOCK_READ = %d\n", BLOCK_READ);
//		#endif
		// 特徴データ全体をdatabaseに読み込まないで，照合のときに読み込む
	} else if(search_mode == 4) {
		// search_mode = 3;
	} else {
		prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);
		if(search_mode == 0) {
			alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
			sort_database2();
		} else if(search_mode == 3) {
			sort_database_self();
		}
	}
	sk_c = sk_c_0 =0;
	clock_gettime(CLOCK_REALTIME, &tp1);

//	for(q = query_from; q < query_from + query_num; q++) {
	for(int qq = 0; qq < query_num; qq++) {
//		fprintf(stderr, "q = %d\n", qq);
		int q = qq + q_num_from;
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE_c2 qu, qu2;
		struct_que_c2_n *que = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n));
		#if ! defined(BLOCK_READ) && ! defined(NUM_THREADS)
		dist_type dist = 0;
		#endif
		int i = 0, j, k;

		ans[q].dist = UINT_MAX;
		clock_gettime(CLOCK_REALTIME, &tp3);
		bound_with_idx(query[q].data, bd_idx);
		if(lb >= 5) {
			p = (double)lb / 10;
			for(i = 0; i < PJT_DIM; i++)
				bd_idx[i].bd = pow(bd_idx[i].bd, p);
		}
		qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

		k = 0;

		#if defined(BLOCK_READ) || defined(NUM_THREADS)
		num_idx = 0;
		#endif
		
		s = query[q].sketch;
		for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
		#if defined(BLOCK_READ) || defined(NUM_THREADS)
			ftr_idx_lst[num_idx++] = j;
			continue;
		#else
			#if defined(SEARCH_ON_SECONDARY) // if(search_mode == 6 || search_mode == 5) 
			if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
				fprintf(stderr, "get_ftr error at %d\n", j);
				exit(0);
			}
			if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
				fprintf(stderr, "cannot read sub_fb at %d\n", j);
				exit(0);
			}
			dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
			#else
			if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 3 || search_mode == 4) {
				dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			#endif
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		#endif
		}

		s = s ^ (1 <<  bd_idx[0].idx);
		for(j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
		#if defined(BLOCK_READ) || defined(NUM_THREADS)
			ftr_idx_lst[num_idx++] = j;
			continue;
		#else
			#if defined(SEARCH_ON_SECONDARY) // if(search_mode == 6 || search_mode == 5) 
			if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
				fprintf(stderr, "get_ftr error at %d\n", j);
				exit(0);
			}
			if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
				fprintf(stderr, "cannot read sub_fb at %d\n", j);
				exit(0);
			}
			dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
			#else
			if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else if(search_mode == 3 || search_mode == 4) {
				dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			#endif
			if(dist < ans[q].dist) {
				ans[q].dist = dist;
				ans[q].idx = idx[j];
			}
		#endif
		}

		make_empty_que_c2_n(que);

		// enq pattern of 0...10
		qu.cursor = new_que_e2_n(que);
		qu.key = bd_idx[1].bd;
		que->details[qu.cursor].sk = query[q].sketch ^ (1 << bd_idx[1].idx);
		que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
		enq_c2_n(&qu, que);		

		while(deq_c2_n(&qu, que) && k < K) {
			s = que->details[qu.cursor].sk;
			sk_c_0 += (bkt[s] == bkt[s + 1]);
			for(j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			#if defined(BLOCK_READ) || defined(NUM_THREADS)
				ftr_idx_lst[num_idx++] = j;
				continue;
			#else
				#if defined(SEARCH_ON_SECONDARY) // if(search_mode == 6 || search_mode == 5) 
				if(get_ftr(dbfh, db_header, j, ftr_data) == 0) {
					fprintf(stderr, "get_ftr error at %d\n", j);
					exit(0);
				}
				if(!READ(dbfh, &data_id, sizeof(sub_fb_type))) {
					fprintf(stderr, "cannot read sub_fb at %d\n", j);
					exit(0);
				}
				dist = DISTANCE_2(ftr_data, query[q].data, db_header->data_dim, ans[q].dist);
				#else
				if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else if(search_mode == 3 || search_mode == 4) {
					dist = DISTANCE_2(database[j], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				#endif
				if(dist < ans[q].dist) {
					ans[q].dist = dist;
					ans[q].idx = idx[j];
				}
			#endif
			}

			switch(que->details[qu.cursor].pt & 15) {
			case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
			case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
				{
					int m = lsb_pos(que->details[qu.cursor].pt);
					if(m > 0 && m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
						// Y010^m -> Y10^{m+1}
						qu2.cursor = new_que_e2_n(que);
						qu2.key = qu.key + bd_idx[m + 1].bd - bd_idx[m].bd;
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[m + 1].idx)) ^ (1 << bd_idx[m].idx);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
						// Y010^m -> Y010^{m-1}1
						qu.key = qu.key + bd_idx[0].bd;
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
						enq_c2_n(&qu2, que);
					} else {
						qu.key = qu.key + bd_idx[0].bd;
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
					}
				}
				break;
			case 4:  // X0100 -> enq(X0101) and enq(X1000)
				// X1000
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + bd_idx[3].bd - bd_idx[2].bd;
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[3].idx)) ^ (1 << bd_idx[2].idx);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
				// X0101
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 1:  // X0001 -> enq(X0010)
			case 5:  // X0101 -> enq(X0110)
			case 9:  // X1001 -> enq(X1010)
			case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
				qu.key = qu.key + bd_idx[1].bd - bd_idx[0].bd;
				que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[1].idx)) ^ (1 << bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 2:  // X0010 -> enq(X0011) and enq(X0100)
			case 10: // X1010 -> enq(X1011) and enq(X1100)
				// X0100 and X1100
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key +  bd_idx[2].bd -  bd_idx[1].bd;
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[2].idx)) ^ (1 << bd_idx[1].idx);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
				// X0011 and X1011
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 6:  // X0110 -> enq(X0111)
			case 12: // X1100 -> enq(X1101)
			case 14: // X1110 -> enq(10111)
				qu.key = qu.key + bd_idx[0].bd;
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 3:  // X0011
			case 7:  // X0111
			case 11: // X1011
			case 15: // X1111 -> nothing to do
				break;
			}
			i++;
		}
		sk_c += i;
		free(que);

		#ifdef NUM_THREADS
		int min_t = 0;
		#endif

	#ifdef WITHOUT_SFTR
			// SFTR無しでの検索は，列挙のみで何もしない．
	#else
		if(search_mode == 6) {
			// search_mode == 6 のときは，列挙のみなので，何もしない
		} else {
			// search_mode == 4 のときは，オンメモリ（ソート済み）データを用いる
		#if ! defined(BLOCK_READ) && ! defined(NUM_THREADS)
			// 何もしない
		#else
			#if defined(BLOCK_READ) && defined(NUM_THREADS)
			for(int t = 0; t < nt; t++) {
				hd_pool[t]->read_in = 0;
				hd_pool[t]->next = 0;
				hd_pool[t]->ans.dist = UINT_MAX;
			}
			#elif defined(BLOCK_READ) && ! defined(NUM_THREADS)
			handle->read_in = 0;
			handle->next = 0;
			handle->ans.dist = UINT_MAX;
			#else
			for(int t = 0; t < nt; t++) {
				temp_ans[t].dist = UINT_MAX;
			}
			#endif

			#ifdef NUM_THREADS
			#pragma omp parallel for
			#endif
			for(k = 0; k < K; k++) {
				ftr_type ftr;
				int ftr_idx;
				dist_type dist;

				#ifdef NUM_THREADS
				int t = omp_get_thread_num(); // 実行中のスレッド番号を取得
				#ifdef BLOCK_READ
				ftr_handle *handle = hd_pool[t];
				#else
				#endif
				#endif

				#ifdef BLOCK_READ
				// ブロック内に読み込んだデータがなくなったら，つぎのブロックをバッファーに読み込む
				if(handle->next >= handle->read_in) {
					get_ftr_block(handle, ftr_idx_lst, k, num_idx);
				}
				// バッファーの次のデータとidxを受け取る
				ftr = handle->ftr_buff[handle->next].ftr;
				ftr_idx = handle->idx[handle->next];
				handle->next++;
				dist = DISTANCE_2(ftr, query[q].data, db_header->data_dim, handle->ans.dist);
				if(dist < handle->ans.dist) {
					handle->ans.dist = dist;
					handle->ans.idx = idx[ftr_idx];
				}
				#else
				ftr_idx = ftr_idx_lst[k];
				ftr = database[ftr_idx];
				dist = DISTANCE_2(ftr, query[q].data, db_header->data_dim, temp_ans[t].dist);
				if(dist < temp_ans[t].dist) {
					temp_ans[t].dist = dist;
					temp_ans[t].idx = idx[ftr_idx];
				}
				#endif
			}
			// 各スレッドが見つけた暫定解から距離最小のものを選ぶ
			#ifdef NUM_THREADS
				dist_type min_dist = UINT_MAX; 
				for(int t = 0; t < nt; t++) {
					#ifdef BLOCK_READ
					if(hd_pool[t]->ans.dist < min_dist) {
						min_dist = hd_pool[t]->ans.dist;
						min_t = t;
					}
					#else
					if(temp_ans[t].dist < min_dist) {
						min_dist = temp_ans[t].dist;
						min_t = t;
					}
					#endif
				}
				#ifdef BLOCK_READ
				ans[q] = hd_pool[min_t]->ans;
				#else
				ans[q] = temp_ans[min_t];
				#endif
			#else
				ans[q] = handle->ans;
			#endif
		#endif // BLOCK_READ || NUM_THREADS
		}
	#endif
		clock_gettime(CLOCK_REALTIME,&tp4);
		sec = tp4.tv_sec - tp3.tv_sec;
		nsec = tp4.tv_nsec - tp3.tv_nsec;
		if(nsec < 0){
			sec--;
			nsec += 1000000000L;
		}
	#ifndef NO_INTERMEDIATE_RESULT
//		#ifdef BLOCK_READ
//			#ifdef NUM_THREADS
//			printf("q = %d, dist = %d, idx = %d, %ld.%09ld (sec), min_t = %d\n", q, ans[q].dist, ans[q].idx, sec, nsec, min_t);
//			#else
//			printf("q = %d, dist = %d, idx = %d, %ld.%09ld (sec)\n", q, ans[q].dist, ans[q].idx, sec, nsec);
//			#endif
//		#else
		printf("q = %d, dist = %d, idx = %d, %ld.%09ld (sec)\n", q, ans[q].dist, ans[q].idx, sec, nsec);
//		#endif
	#endif
	}
	#ifdef BLOCK_READ
	#ifndef NUM_THREADS
		free_ftr_handle(handle);
	#else
		for(int t = 0; t < nt; t++) {
			free_ftr_handle(hd_pool[t]);
		}
	#endif
	#endif
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
}

void L1_search_16bit_2_kNN(int search_mode) { // 列挙による高速化あり
	double p;
	int i, q, k;
	int *idx = NULL;
	int *bkt = NULL;
	int b_idx[PJT_DIM];
	dist_type dist = 0;
	sketch_type s;
	QUE qu, qu2;
	int skipped = 0, total_skipped = 0;

	prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);

	if(database2 == NULL) {
		alloc_database2(db_header->element_size, db_header->data_dim, db_header->data_num);
	}
	for(i = 0; i < db_header->data_num; i++) {
		memcpy(database2[i], database[idx[i]], FTR_DIM);
	}
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;

// kNN search のための書き換え部分（ここから）
	kNN_answer = (kNN_buffer **)malloc(sizeof(kNN_buffer *) * query_num);
	for(i = 0; i < query_num; i++) {
		kNN_answer[i] = new_kNN_buffer();
	}
// kNN search のための書き換え部分（ここまで）

	sk_c = sk_c_0 =0;
	
	init_que();
	
	clock_gettime(CLOCK_REALTIME, &tp1);
	
	for(q = 0; q < query_num; q++) {
		bound(query[q].data, bd);

		if(lb >= 5) {
			p = (double)lb / 10;
			for(i = 0; i < PJT_DIM; i++)
				bd[i] = pow(bd[i], p);
		}

		for(i = 0; i < PJT_DIM; i++)
			b_idx[i] = i;
		qsort(b_idx, PJT_DIM, sizeof(int), comp_b);
		s = query[q].sketch;
		k = 0;
		skipped = 0;
		for(int j = bkt[s]; j < bkt[s + 1] && k < K; j++, k++) {
			if(search_mode == 2) {
				// nothing to do
			} else if(search_mode == 1) {
				dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
			} else {
				dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
			}
			ans[q].dist = dist;
			ans[q].idx = idx[j];
			dist = push_kNN_buffer(&ans[q], kNN_answer[q]);
			ans[q].dist = dist;
		}

		qu.key = bd[b_idx[0]];
		qu.sk = s ^ (1 << b_idx[0]);
		qu.idx = 1;
		enq(&qu);
		i = 1;

		while(deq(&qu) && k < K) {
			sk_c_0 += (bkt[qu.sk] == bkt[qu.sk + 1]);
			for(int j = bkt[qu.sk]; j < bkt[qu.sk + 1] && k < K; j++, k++) {
				if(search_mode == 2) {
					// nothing to do
				} else if(search_mode == 1) {
					dist = DISTANCE_2(database[idx[j]], query[q].data, db_header->data_dim, ans[q].dist);
				} else {
					dist = DISTANCE_2(database2[j], query[q].data, db_header->data_dim, ans[q].dist);
				}
				ans[q].dist = dist;
				ans[q].idx = idx[j];
				dist = push_kNN_buffer(&ans[q], kNN_answer[q]);
				ans[q].dist = dist;
			}
			
			if(qu.idx < PJT_DIM) {
				s = qu.sk ^ (1 << b_idx[qu.idx]);
				qu2.key = qu.key + bd[b_idx[qu.idx]];
				qu2.sk = s;
				qu2.idx = qu.idx + 1;
#ifdef SKIP_ENQ
				if(SKIP_CONDITION(qu2.key, ans[q].dist)) {
					skipped++;
				} else {
					enq(&qu2);
				}
#else
				enq(&qu2);
#endif
				qu2.key = qu.key + bd[b_idx[qu.idx]] - bd[b_idx[qu.idx - 1]];
				qu2.sk = s ^ (1 << b_idx[qu.idx - 1]);
				qu2.idx = qu.idx + 1;
#ifdef SKIP_ENQ
				if(SKIP_CONDITION(qu2.key, ans[q].dist)) {
					skipped++;
				} else {
					enq(&qu2);
				}
#else
				enq(&qu2);
#endif
			}
			i++;
		}
		sk_c += i;
		qsize = 0; // priority queueの初期化
		total_skipped += skipped;
		flush_kNN_buffer(kNN_answer[q]);
	}
	
	clock_gettime(CLOCK_REALTIME,&tp2);
	fprintf(stderr, "\n");
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	fprintf(stdout, "%ld.%09ld\n", sec, nsec);
	fprintf(stderr,"total_skipped = %d\n", total_skipped);
}

int countBit16(sketch_type n) {
	n = (n & 0x5555) + ((n >> 1) & 0x5555);
	n = (n & 0x3333) + ((n >> 2) & 0x3333);
	n = (n + (n >> 4)) & 0x0F0F;
	return (n + (n >> 8)) & 0xFF;
}

int countBit32(sketch_type n) {
	n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	n = (n + (n >> 4)) & 0x0F0F0F0F;
	n += n >> 8;
	return (n + (n >> 16)) & 0xFF;
}

void search_K_multi(int multi_K[]) {
	int i, q;

	prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;

//	omp_set_num_threads(4);
#ifdef TEST_PARALLEL
	#pragma omp parallel for
#endif
	for(q = 0; q < query_num; q++) {
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE qu, qu2;
		struct_que *que;
		int i, k = 0;

		que = (struct_que *)malloc(sizeof(struct_que));
		if(lb == 1 || lb >= 5) { // L1
			bound_with_idx(query[q].data, bd_idx);
			if(lb >= 5) {
				p = (double)lb / 10;
				for(i = 0; i < PJT_DIM; i++)
					bd_idx[i].bd = pow(bd_idx[i].bd, p);
			}
			qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

			s = query[q].sketch;
			for(int j = bkt[s], k = 0; j < bkt[s + 1] && k < db_header->data_num; j++, k++) {
				if(idx[j] == query[q].nearest) {
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					k = db_header->data_num; // 検索終了
				}
			}

			que->qsize = 0;
			qu.key = bd_idx[0].bd;
			qu.sk = s ^ (1 << bd_idx[0].idx);
			qu.idx = 1;
			ENQ_P(&qu, que);

			while(DEQ_P(&qu, que) && k < db_header->data_num) {
				for(int j = bkt[s]; j < bkt[qu.sk + 1] && k < db_header->data_num; j++, k++) {
					if(idx[j] == query[q].nearest) {
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						k = db_header->data_num; // 検索終了
					}
				}
				if(qu.idx < PJT_DIM) {
					s = qu.sk ^ (1 << bd_idx[qu.idx].idx);
					qu2.key = qu.key + bd_idx[qu.idx].bd;
					qu2.sk = s;
					qu2.idx = qu.idx + 1;
					ENQ_P(&qu2, que);

					qu2.key = qu.key + bd_idx[qu.idx].bd - bd_idx[qu.idx - 1].bd;
					qu2.sk = s ^ (1 << bd_idx[qu.idx - 1].idx);
					qu2.idx = qu.idx + 1;
					ENQ_P(&qu2, que);
				}
			}
		}
		free(que);
	}

	qsort(ans, query_num, sizeof(answer), comp_ans);
	for(i = 0; i < 1000; i++) {
		multi_K[i] = (int)(ans[(int)(query_num * i / 1000)].dist);
	}
}

void search_K_multi_c2_n(int multi_K[]) {
	int i, q;

	prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);

#ifdef TEST_PARALLEL
	int nt = omp_get_max_threads();
	#pragma omp parallel for
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;

	static struct_que_c2_n *que_pool = NULL;
	int used_q[nt];
	if(que_pool == NULL) {
		fprintf(stderr, "Malloc for que requests (%ld * %d = %ld is required).\n", sizeof(struct_que_c2_n), nt, sizeof(struct_que_c2_n) * nt);
		que_pool = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n) * nt);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed\n");
			exit(0);
		}
		for(q = 0; q < nt; q++) {
			make_empty_que_c2_n(&que_pool[q]);
			used_q[q] = 0;
		}
	}
#else
	struct_que_c2_n *que = NULL;
	que = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n));
	if(que == NULL) {
		fprintf(stderr, "Malloc for que failed\n");
		exit(0);
	}
	make_empty_que_c2_n(que);
#endif

#ifdef TEST_PARALLEL
	#pragma omp parallel for
#endif
	for(q = 0; q < query_num; q++) {
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE_c2 qu, qu2;
		#ifdef TEST_PARALLEL
		struct_que_c2_n *que = NULL;
		#endif
		int i, k = 0;

		if(lb == 1 || lb >= 5) { // L1
			bound_with_idx(query[q].data, bd_idx);
			if(lb >= 5) {
				p = (double)lb / 10;
				for(i = 0; i < PJT_DIM; i++)
					bd_idx[i].bd = pow(bd_idx[i].bd, p);
			}
			qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

			s = query[q].sketch;
			for(int j = bkt[s], k = 0; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				}
			}

			s = s ^ (1 << bd_idx[0].idx);
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				}
			}

			#ifdef TEST_PARALLEL
			int tn;
			#pragma omp critical 
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
			#endif

			make_empty_que_c2_n(que);

			// enq pattern of 0...10
			qu.cursor = new_que_e2_n(que);
			qu.key = bd_idx[1].bd;
			que->details[qu.cursor].sk = query[q].sketch ^ (1 << bd_idx[1].idx);
			que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
			enq_c2_n(&qu, que);

			while(deq_c2_n(&qu, que)) {
				s = que->details[qu.cursor].sk;
				for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
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
							qu2.key = qu.key + bd_idx[m + 1].bd - bd_idx[m].bd;
							que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[m + 1].idx)) ^ (1 << bd_idx[m].idx);
							que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
							// Y010^m -> Y010^{m-1}1
							qu.key = qu.key + bd_idx[0].bd;
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							enq_c2_n(&qu, que);
							enq_c2_n(&qu2, que);
						} else {
							qu.key = qu.key + bd_idx[0].bd;
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							enq_c2_n(&qu, que);
						}
					}
					break;
				case 4:  // X0100 -> enq(X0101) and enq(X1000)
					// X1000
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key +  bd_idx[3].bd -  bd_idx[2].bd;
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[3].idx)) ^ (1 << bd_idx[2].idx);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
					// X0101
					qu.key = qu.key + bd_idx[0].bd;
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
					break;
				case 1:  // X0001 -> enq(X0010)
				case 5:  // X0101 -> enq(X0110)
				case 9:  // X1001 -> enq(X1010)
				case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
					qu.key = qu.key + bd_idx[1].bd - bd_idx[0].bd;
					que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[1].idx)) ^ (1 << bd_idx[0].idx);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					break;
				case 2:  // X0010 -> enq(X0011) and enq(X0100)
				case 10: // X1010 -> enq(X1011) and enq(X1100)
					// X0100 and X1100
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key + bd_idx[2].bd - bd_idx[1].bd;
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[2].idx)) ^ (1 << bd_idx[1].idx);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
					// X0011 and X1011
					qu.key = qu.key + bd_idx[0].bd;
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
					break;
				case 6:  // X0110 -> enq(X0111)
				case 12: // X1100 -> enq(X1101)
				case 14: // X1110 -> enq(10111)
					qu.key = qu.key + bd_idx[0].bd;
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
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
			#ifdef TEST_PARALLEL
			#pragma omp critical 
			{
				used_q[tn] = 0;
			}
			#endif
END_SEARCH:	;
		}
	}

	qsort(ans, query_num, sizeof(answer), comp_ans);
	for(i = 0; i < 1000; i++) {
		multi_K[i] = (int)(ans[(int)(query_num * i / 1000)].dist);
	}
}

void search_K_multi_c2_n_new_bucket(int multi_K[]) {
	int i, q;

//	prepare_bucket_db(1, sketch, db_header->data_num, &idx, &bkt);

#ifdef TEST_PARALLEL
	int nt = omp_get_max_threads();
	if(PJT_DIM > 20) {
		nt /= 4;
		omp_set_num_threads(nt);
	}
#endif
	#pragma omp parallel for
	for(i = 0; i < query_num; i++)
		ans[i].dist = UINT_MAX;

#ifdef TEST_PARALLEL
	static struct_que_c2_n *que_pool = NULL;
//	int used_q[nt];
	if(que_pool == NULL) {
		fprintf(stderr, "Malloc for que requests (%ld * %d = %ld is required).\n", sizeof(struct_que_c2_n), nt, sizeof(struct_que_c2_n) * nt);
		que_pool = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n) * nt);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed\n");
			exit(0);
		}
//		for(q = 0; q < nt; q++) {
//			make_empty_que_c2_n(&que_pool[q]);
//			used_q[q] = 0;
//		}
	}
#else
	struct_que_c2_n *que = (struct_que_c2_n *)malloc(sizeof(struct_que_c2_n));

#endif

#ifdef TEST_PARALLEL
	printf("parallel: nt = %d\n", nt);
	#pragma omp parallel for
#endif
	for(q = 0; q < query_num; q++) {
		sketch_type s;
		double p;
		bd_type bd_idx[PJT_DIM];
		QUE_c2 qu, qu2;
		int qn = q + q_num_from;
#ifdef TEST_PARALLEL
		struct_que_c2_n *que = &que_pool[omp_get_thread_num()];
#endif
		int i, j, k;

		if(lb == 1 || lb >= 5) { // L1
			bound_with_idx(query[qn].data, bd_idx);
			if(lb >= 5) {
				p = (double)lb / 10;
				for(i = 0; i < PJT_DIM; i++)
					bd_idx[i].bd = pow(bd_idx[i].bd, p);
			}
			qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

			s = query[qn].sketch;
			for(k = 0, j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[qn].nearest) {
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				}
			}

			s = s ^ (1 << bd_idx[0].idx);
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[qn].nearest) {
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				}
			}

/*
			#ifdef TEST_PARALLEL
			int tn = omp_get_thread_num();
			que = &que_pool[tn];
			
			int tn;
			#pragma omp critical 
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

			#endif
*/			
			make_empty_que_c2_n(que);

			// enq pattern of 0...10
			qu.cursor = new_que_e2_n(que);
			qu.key = bd_idx[1].bd;
			que->details[qu.cursor].sk = query[qn].sketch ^ (1 << bd_idx[1].idx);
			que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
			enq_c2_n(&qu, que);

			while(deq_c2_n(&qu, que)) {
				s = que->details[qu.cursor].sk;
				for(j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[qn].nearest) {
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
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
							qu2.key = qu.key + bd_idx[m + 1].bd - bd_idx[m].bd;
							que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[m + 1].idx)) ^ (1 << bd_idx[m].idx);
							que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
							// Y010^m -> Y010^{m-1}1
							qu.key = qu.key + bd_idx[0].bd;
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							enq_c2_n(&qu, que);
							enq_c2_n(&qu2, que);
						} else {
							qu.key = qu.key + bd_idx[0].bd;
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							enq_c2_n(&qu, que);
						}
					}
					break;
				case 4:  // X0100 -> enq(X0101) and enq(X1000)
					// X1000
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key +  bd_idx[3].bd -  bd_idx[2].bd;
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[3].idx)) ^ (1 << bd_idx[2].idx);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
					// X0101
					qu.key = qu.key + bd_idx[0].bd;
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0].idx);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
					break;
				case 1:  // X0001 -> enq(X0010)
				case 5:  // X0101 -> enq(X0110)
				case 9:  // X1001 -> enq(X1010)
				case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
					qu.key = qu.key + bd_idx[1].bd - bd_idx[0].bd;
					que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[1].idx)) ^ (1 << bd_idx[0].idx);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					break;
				case 2:  // X0010 -> enq(X0011) and enq(X0100)
				case 10: // X1010 -> enq(X1011) and enq(X1100)
					// X0100 and X1100
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key + bd_idx[2].bd - bd_idx[1].bd;
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[2].idx)) ^ (1 << bd_idx[1].idx);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
					// X0011 and X1011
					qu.key = qu.key + bd_idx[0].bd;
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
					break;
				case 6:  // X0110 -> enq(X0111)
				case 12: // X1100 -> enq(X1101)
				case 14: // X1110 -> enq(10111)
					qu.key = qu.key + bd_idx[0].bd;
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0].idx);
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
/*
			#ifdef TEST_PARALLEL
			#pragma omp critical 
			{
				used_q[tn] = 0;
			}
			#endif
*/
END_SEARCH:	;
		}
	}

	qsort(ans, query_num, sizeof(answer), comp_ans);

	for(i = 0; i < 1000; i++) {
		multi_K[i] = (int)(ans[(int)(query_num * i / 1000)].dist);
	}
}

static int bucket_database_flag = 0; // 0 = not prepared, 1 = fixed, 2 = variable

int prepare_bucket_db(int flag, sketch_type *sk, int data_num, int **idx, int **bkt)
// version 2.0 (using arrays only idx and bucket)
{
	int i;
	int *bucket = NULL;

	if(bucket_database_flag == 0) {
		if(*idx == NULL) *idx = (int *)malloc(sizeof(int) * data_num);
		if(*bkt == NULL) *bkt = (int *)calloc(BIT + 2, sizeof(int));
	} else if(bucket_database_flag == 1) {
		fprintf(stderr, "Bucket database (ver. 2) is already prepared.\n");
		return flag;
	}

	bucket = *bkt + 1;
	for(i = 0; i < BIT + 1; i++) bucket[i] = 0;
	for(i = 0; i < data_num; i++) bucket[sk[i] + 1]++;
	for(i = 0; i < BIT; i++) bucket[i + 1] += bucket[i];
	for(i = 0; i < data_num; i++) (*idx)[bucket[sk[i]]++] = i;
	(*bkt)[0] = 0;

	bucket_database_flag = flag;

	return flag;
}

void write_bucket_db(char *filename, int data_num, int *idx, int *bkt)
{
	FILE *fp;
	int num_elements, num_nonempty_buckets = 0;
	sketch_type s;

	if((fp = fopen(filename, "wb"))  == NULL) {
		fprintf(stderr, "Write open bucket file error, file name = %s\n", filename);
		exit(0);
	}
	if(fwrite(&data_num, sizeof(int), 1, fp) != 1) {  // data_num を書き出す
		fprintf(stderr, "fwrite error (data_num) file = %s\n", filename);
		fclose(fp);
		exit(0);
	}
	if(fwrite(idx, sizeof(int) * db_header->data_num, 1, fp) != 1) {  // idx[data_num] を書き出す
		fprintf(stderr, "fwrite error (idx) file = %s\n", filename);
		fclose(fp);
		exit(0);
	}
	for(s = 0; s < BIT; s++) {
		num_nonempty_buckets += (bkt[s + 1] - bkt[s] > 0);
	}
	if(fwrite(&num_nonempty_buckets, sizeof(int), 1, fp) != 1) {  // num_nonempty_buckets を書き出す
		fprintf(stderr, "fwrite error (num_nonempty_buckets) file = %s\n", filename);
		fclose(fp);
		exit(0);
	}
	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", num_nonempty_buckets, (double)data_num / num_nonempty_buckets);
	for(s = 0; s < BIT; s++) {
		num_elements = bkt[s + 1] - bkt[s];
		if(num_elements > 0) {
			if(fwrite(&s, sizeof(sketch_type), 1, fp) != 1) {  // sketch s を書き出す
				fprintf(stderr, "fwrite error (sketch) file = %s\n", filename);
				exit(0);
				return;
			}
			if(fwrite(&num_elements, sizeof(int), 1, fp) != 1) {  // num_elements を書き出す
				fprintf(stderr, "fwrite error (num_elements) file = %s\n", filename);
				fclose(fp);
				exit(0);
			}
		}
	}
	fclose(fp);
	return;
}

int read_bucket_db(char *filename, int data_num, int **idx, int **bkt)
{
	FILE *fp;
	int num;
	int num_elements, num_nonempty_buckets = 0;
	sketch_type s, s_next;
	
	if((fp = fopen(filename, "rb"))  == NULL) {
		fprintf(stderr, "Read open bucket file error, file name = %s\n", filename);
		return 0;
	}

	if(fread(&num, sizeof(int), 1, fp) != 1) {  // ファイルに書かれている data_num を読み込む
		fprintf(stderr, "fread error (data_num) file = %s\n", filename);
		fclose(fp);
		return 0;
	}
	if(data_num != 0 && data_num != num) { // データ数についてのチェック（ただし，data_num == 0 のときは，ノーチェック）
		fprintf(stderr, "data_num confliction (request = %d, in_file = %d)\n", data_num, num);
		fclose(fp);
		return 0;
	}
	data_num = num;
	
	if(*idx == NULL) *idx = (int *)malloc(sizeof(int) * data_num);
	if(*bkt == NULL) *bkt = (int *)calloc(BIT + 2, sizeof(int));
	
	if(fread(*idx, sizeof(int) * data_num, 1, fp) != 1) {  // idx[data_num] を読み込む
		fprintf(stderr, "fread error (idx, size = %ld) file = %s\n", sizeof(int) * db_header->data_num, filename);
		fclose(fp);
		return 0;
	}
	if(fread(&num_nonempty_buckets, sizeof(int), 1, fp) != 1) {  // ファイルに書かれている num_nonempty_buckets を読み込む
		fprintf(stderr, "fread error (data_num) file = %s\n", filename);
		fclose(fp);
		return 0;
	}
	for(s = 0; s < BIT + 2; s++) (*bkt)[s] = 0;
	num = s_next = 0;
	for(int i = 0; i < num_nonempty_buckets; i++) {
		if(fread(&s, sizeof(sketch_type), 1, fp) != 1) {  // sketch を読み込む
			fprintf(stderr, "fread error (sketch) file = %s\n", filename);
			fclose(fp);
			return 0;
		}
		if(fread(&num_elements, sizeof(int), 1, fp) != 1) {  // num_elements を読み込む
			fprintf(stderr, "fread error (sketch) file = %s\n", filename);
			fclose(fp);
			return 0;
		}
		while(s_next <= s) {
			(*bkt)[s_next++] = num;
		}
		num += num_elements;
	}
	while(s_next < BIT + 2) {
		(*bkt)[s_next++] = num;
	}
	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", num_nonempty_buckets, (double)data_num / num_nonempty_buckets);
	fclose(fp);
	return num;
}

double precision_resample(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	// Prool of priority queue for parallel
	static struct_que *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que *)malloc(sizeof(struct_que) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		QUE qu, qu2;
		struct_que *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均

		int tn;
		#pragma omp critical 
		{
			for(tn = 0; tn < 16; tn++) {
				if(used_q[tn] == 0) {
					que = &que_pool[tn];
					used_q[tn] = 1;
					break;
				}
				if(tn == 16) {
					fprintf(stderr, "No free que in que_pool\n");
					exit(0);
				}
			}
		}

		if(lb == 1 || lb >= 5) { // L1
			s = query[q].sketch;
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
//			if(debug_flag) printf("q, %d, s, %d, k, %d, min_k, %d\n", q, s, k, min_k);
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
//					sum += log(SCORE_FACTOR * db_header->data_num / (min_k + 1));
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto FREE_Q; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			que->qsize = 0;
			qu.key = query[q].bd[query[q].idx[0]];
			qu.sk = s ^ (1 << query[q].idx[0]);
			qu.idx = 1;
			ENQ_P(&qu, que);

			while(DEQ_P(&qu, que)) {
				if(qu.key > score) {
					if(k >= trunc_K) goto FREE_Q;
					score = qu.key;
//					min_k = k;
				}
//				if(debug_flag) printf("q, %d, s, %d, k, %d, min_k, %d\n", q, qu.sk, k, min_k);
				mid_k = k +  (bkt[qu.sk + 1] - bkt[qu.sk]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[qu.sk]; j < bkt[qu.sk + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
//						sum += log(SCORE_FACTOR * db_header->data_num / (min_k + 1));
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[qu.sk + 1] - j;
						break;
					}
				}
				if(qu.idx < PJT_DIM) {
					s = qu.sk ^ (1 << query[q].idx[qu.idx]);
					qu2.key = qu.key + query[q].bd[query[q].idx[qu.idx]];
					qu2.sk = s;
					qu2.idx = qu.idx + 1;
					ENQ_P(&qu2, que);

					qu2.key = qu.key + query[q].bd[query[q].idx[qu.idx]] - query[q].bd[query[q].idx[qu.idx - 1]];
					qu2.sk = s ^ (1 << query[q].idx[qu.idx - 1]);
					qu2.idx = qu.idx + 1;
					ENQ_P(&qu2, que);
				}
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
			}
		}
	}
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// 再サンプリング(大きさ = resampling_size)の質問に対する精度（delayed enq）
double precision_resample_delayed(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	// Prool of priority queue for parallel
	static struct_que *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que *)malloc(sizeof(struct_que) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		QUE qu, qu2;
		struct_que *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif

			s = query[q].sketch;
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
				if(k >= trunc_K) goto END_SEARCH;
			}
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
			{
				for(tn = 0; tn < 16; tn++) {
					if(used_q[tn] == 0) {
						que = &que_pool[tn];
						used_q[tn] = 1;
						break;
					}
					if(tn == 16) {
						fprintf(stderr, "No free que in que_pool\n");
						exit(0);
					}
				}
			}

			que->qsize = 0;
			qu.key = query[q].bd[query[q].idx[1]];
			qu.sk = 1 << query[q].idx[1];
			qu.idx = 2;
			ENQ_P(&qu, que);

			while(DEQ_P(&qu, que)) {
				s = query[q].sketch ^ qu.sk;
				if(qu.key > score) {
					if(k >= trunc_K) goto FREE_Q;
					score = qu.key;
				}
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; // k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}
				if(qu.idx < PJT_DIM) {
					qu2.key = qu.key + query[q].bd[query[q].idx[qu.idx]] - query[q].bd[query[q].idx[qu.idx - 1]];
					qu2.sk = (qu.sk ^ (1 << query[q].idx[qu.idx - 1])) ^ (1 << query[q].idx[qu.idx]);
					qu2.idx = qu.idx + 1;
					ENQ_P(&qu2, que);

					if(!(qu.sk & (1 << query[q].idx[qu.idx - 2]))) {
						qu.key = qu.key + query[q].bd[query[q].idx[qu.idx - 2]];
						qu.sk = qu.sk ^ (1 << query[q].idx[qu.idx - 2]);
						ENQ_P(&qu, que);
					}
				}
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
			}
END_SEARCH: ;
		}
	}
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// カーソルを使った優先度付き待ち行列を使用するバージョン（delayed enq を用いている）
double precision_resample_cursor(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	static struct_que_c *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que_c *)malloc(sizeof(struct_que_c) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
		for(q = 0; q < 16; q++) {
			make_empty_que_c(&que_pool[q]);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		int qu, qu2;
		struct_que_c *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
			s = query[q].sketch;
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j; 
					break;
				}
			}
			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
				if(k >= trunc_K) goto END_SEARCH;
			}
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
			{
				for(tn = 0; tn < 16; tn++) {
					if(used_q[tn] == 0) {
						que = &que_pool[tn];
						used_q[tn] = 1;
						break;
					}
					if(tn == 16) {
						fprintf(stderr, "No free que in que_pool\n");
						exit(0);
					}
				}
			}

			make_empty_que_c(que);

			qu = new_que_e(que);
			que->que_e[qu].key = query[q].bd[query[q].idx[1]];
			que->que_e[qu].sk = 1 << query[q].idx[1];
			que->que_e[qu].idx = 2;
			enq_c(qu, que);

			while(deq_c(&qu, que)) {
				s = query[q].sketch ^ que->que_e[qu].sk;
				if(que->que_e[qu].key > score) {
					if(k >= trunc_K) goto FREE_Q;
					score = que->que_e[qu].key;
				}
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}
				if(que->que_e[qu].idx < PJT_DIM) {
					qu2 = new_que_e(que);
					que->que_e[qu2].key = que->que_e[qu].key + query[q].bd[query[q].idx[que->que_e[qu].idx]] - query[q].bd[query[q].idx[que->que_e[qu].idx - 1]];
					que->que_e[qu2].sk = (que->que_e[qu].sk ^ (1 << query[q].idx[que->que_e[qu].idx - 1])) ^ (1 << query[q].idx[que->que_e[qu].idx]);
					que->que_e[qu2].idx = que->que_e[qu].idx + 1;
					enq_c(qu2, que);
					if(!(que->que_e[qu].sk & (1 << query[q].idx[que->que_e[qu].idx - 2]))) {
						que->que_e[qu].key = que->que_e[qu].key + query[q].bd[query[q].idx[que->que_e[qu].idx - 2]];
						que->que_e[qu].sk = que->que_e[qu].sk ^ (1 << query[q].idx[que->que_e[qu].idx - 2]);
						enq_c(qu, que);
					}
				}
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
			}
END_SEARCH: ;
		}
	}
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// カーソルを使った優先度付き待ち行列を使用するimamuraバージョン（delayed enq を用いている）
double precision_resample_cursor_2(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif
//	int num_enumerated = 0, num_remained = 0;

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	static struct_que_c2 *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que_c2 *)malloc(sizeof(struct_que_c2) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
		for(q = 0; q < 16; q++) {
			make_empty_que_c2(&que_pool[q]);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
//		printf("q = %d\n", q);
		// local variables for parallel
		sketch_type s;
		QUE_c2 qu, qu2;
		struct_que_c2 *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均
//		int enumerated = 0;

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
			
			s = query[q].sketch;
//			if(debug_flag) printf("q, %d, s, %d, k, %d, min_k, %d\n", q, s, k, min_k);

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j; 
					break;
				}
			}
			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
				if(k >= trunc_K) goto END_SEARCH;
			}

//			if(debug_flag) printf("q, %d, s, %d, k, %d, min_k, %d\n", q, s, k, min_k);
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
			{
				for(tn = 0; tn < 16; tn++) {
					if(used_q[tn] == 0) {
						que = &que_pool[tn];
						used_q[tn] = 1;
						break;
					}
					if(tn == 16) {
						fprintf(stderr, "No free que in que_pool\n");
						exit(0);
					}
				}
			}

			make_empty_que_c2(que);

			// enq pattern of 0...10
			qu.cursor = new_que_e2(que);
			qu.key = query[q].bd[query[q].idx[1]];
			que->details[qu.cursor].sk = 1 << query[q].idx[1];
			que->details[qu.cursor].idx = 2;
			enq_c2(&qu, que);

			while(deq_c2(&qu, que)) {
//START_SEARCH:
//				enumerated++;
				s = query[q].sketch ^ que->details[qu.cursor].sk;
				if(qu.key > score) {
					if(k >= trunc_K) goto FREE_Q;
					score = qu.key;
				}
//				if(q==7) printf("q = %4d, cursor = %6d, key = %6d, k = %6d, min_k = %6d, num = %4d\n", q, qu.cursor, qu.key, k, min_k, sk_num[s]);// getchar();

//				if(debug_flag) printf("q, %d, s, %d, k, %d, min_k, %d\n", q, s, k, min_k);
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}

				if(!(que->details[qu.cursor].sk & (1 << query[q].idx[que->details[qu.cursor].idx - 2]))) {
					if(que->details[qu.cursor].idx < PJT_DIM) {
						qu2.cursor = new_que_e2(que);
						qu2.key = qu.key + query[q].bd[query[q].idx[que->details[qu.cursor].idx]] - query[q].bd[query[q].idx[que->details[qu.cursor].idx - 1]];
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query[q].idx[que->details[qu.cursor].idx - 1])) ^ (1 << query[q].idx[que->details[qu.cursor].idx]);
						que->details[qu2.cursor].idx = que->details[qu.cursor].idx + 1;
						qu.key = qu.key + query[q].bd[query[q].idx[que->details[qu.cursor].idx - 2]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[que->details[qu.cursor].idx - 2]);
						#ifdef USE_TOP
						enq_c2_after_deq(&qu, que);
						#else
						enq_c2(&qu, que);
						#endif
						enq_c2(&qu2, que);
					} else {
						qu.key = qu.key + query[q].bd[query[q].idx[que->details[qu.cursor].idx - 2]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[que->details[qu.cursor].idx - 2]);
						#ifdef USE_TOP
						enq_c2_after_deq(&qu, que);
						#else
						enq_c2(&qu, que);
						#endif
					} 
				} else if(que->details[qu.cursor].idx < PJT_DIM) {
					qu.key = qu.key + query[q].bd[query[q].idx[que->details[qu.cursor].idx]] - query[q].bd[query[q].idx[que->details[qu.cursor].idx - 1]];
					que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query[q].idx[que->details[qu.cursor].idx - 1])) ^ (1 << query[q].idx[que->details[qu.cursor].idx]);
					que->details[qu.cursor].idx = que->details[qu.cursor].idx + 1;
					#ifdef USE_TOP
					enq_c2_after_deq(&qu, que);
					#else
					enq_c2(&qu, que);
					#endif
				} else {
					#ifdef USE_TOP
					deq_c2_del(que);
					#endif
				}
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
//				num_enumerated += enumerated;
//				num_remained += que->qsize;
			}
END_SEARCH: ;
		}
	}
//	printf("enumerated = %d, remained = %d, shurinked = %6.2lf%%\n", num_enumerated, num_remained, (double)(num_enumerated - num_remained) / num_enumerated); fflush(stdout);
//	exit(0);
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// カーソルを使った優先度付き待ち行列を使用するimamuraバージョン（delayed enq を用いている）
// 一部は，樋口のideaによる
double precision_resample_cursor_2_n(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif
//	int num_enumerated = 0, num_remained = 0;

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	int nt = omp_get_max_threads();
	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

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

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
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

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
			
			s = query[q].sketch;

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j; 
					break;
				}
			}
			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
				if(k >= trunc_K) goto END_SEARCH;
			}

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
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
			qu.key = query[q].bd[query[q].idx[1]];
			que->details[qu.cursor].sk = query[q].sketch ^ (1 << query[q].idx[1]);
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
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}
				#ifdef SWITCH8
				switch(que->details[qu.cursor].pt & 15) {
				case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
				case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
					{
						int m = lsb_pos(que->details[qu.cursor].pt);
						if(m > 0 && m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
							// Y010^m -> Y10^{m+1}
							qu2.cursor = new_que_e2_n(que);
							qu2.key = qu.key + query[q].bd[query[q].idx[m + 1]] - query[q].bd[query[q].idx[m]];
							que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query[q].idx[m + 1])) ^ (1 << query[q].idx[m]);
							que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
							// Y010^m -> Y010^{m-1}1
							qu.key = qu.key + query[q].bd[query[q].idx[0]];
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[0]);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							enq_c2_n(&qu, que);
							enq_c2_n(&qu2, que);
						} else {
							qu.key = qu.key + query[q].bd[query[q].idx[0]];
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[0]);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							enq_c2_n(&qu, que);
						}
					}
					break;
				case 4:  // X0100 -> enq(X0101) and enq(X1000)
					// X1000
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key + query[q].bd[query[q].idx[3]] - query[q].bd[query[q].idx[2]];
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query[q].idx[3])) ^ (1 << query[q].idx[2]);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
					// X0101
					qu.key = qu.key + query[q].bd[query[q].idx[0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
					break;
				case 1:  // X0001 -> enq(X0010)
				case 5:  // X0101 -> enq(X0110)
				case 9:  // X1001 -> enq(X1010)
				case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
					qu.key = qu.key + query[q].bd[query[q].idx[1]] - query[q].bd[query[q].idx[0]];
					que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query[q].idx[1])) ^ (1 << query[q].idx[0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					break;
				case 2:  // X0010 -> enq(X0011) and enq(X0100)
				case 10: // X1010 -> enq(X1011) and enq(X1100)
					// X0100 and X1100
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key + query[q].bd[query[q].idx[2]] - query[q].bd[query[q].idx[1]];
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << query[q].idx[2])) ^ (1 << query[q].idx[1]);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
					// X0011 and X1011
					qu.key = qu.key + query[q].bd[query[q].idx[0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
					break;
				case 6:  // X0110 -> enq(X0111)
				case 12: // X1100 -> enq(X1101)
				case 14: // X1110 -> enq(10111)
					qu.key = qu.key + query[q].bd[query[q].idx[0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << query[q].idx[0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					break;
				case 3:  // X0011
				case 7:  // X0111
				case 11: // X1011
				case 15: // X1111 -> nothing to do
					break;
				}
				#else
				switch(que->details[qu.cursor].pt & 3) {
				case 0: // X00 -> enq(X01) and enq(Y10^{m+1}) if X00 = Y010^m
					{
						int m = lsb_pos(que->details[qu.cursor].pt);
						if(m > 0 && m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & MASK(m + 1))) {
							// Y010^m -> Y10^{m+1}
							qu2.cursor = new_que_e2_n(que);
							qu2.key = qu.key + query[q].bd[query[q].idx[m + 1]] - query[q].bd[query[q].idx[m]];
							que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ MASK(query[q].idx[m + 1])) ^ MASK(query[q].idx[m]);
							que->details[qu2.cursor].pt = que->details[qu.cursor].pt + MASK(m);
							// Y010^m -> Y010^{m-1}1
							qu.key = qu.key + query[q].bd[query[q].idx[0]];
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ MASK(query[q].idx[0]);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							#ifdef USE_TOP
							enq_c2_n_after_deq(&qu2, que);
							#else
							enq_c2_n(&qu2, que);
							#endif
							enq_c2_n(&qu, que);
						} else {
							qu.key = qu.key + query[q].bd[query[q].idx[0]];
							que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ MASK(query[q].idx[0]);
							que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
							#ifdef USE_TOP
							enq_c2_n_after_deq(&qu, que);
							#else
							enq_c2_n(&qu, que);
							#endif
						}
					}
					break;
				case 1: // X01 -> enq(X10) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
					qu.key = qu.key + query[q].bd[query[q].idx[1]] - query[q].bd[query[q].idx[0]];
					que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ MASK(query[q].idx[1])) ^ MASK(query[q].idx[0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					#ifdef USE_TOP
					enq_c2_n_after_deq(&qu, que);
					#else
					enq_c2_n(&qu, que);
					#endif
					break;
				case 2: // X10 -> enq(X11) and (enq(x100) only if X10 = X'010)
					if(!(que->details[qu.cursor].pt & 4)) { // X'010 -> enq(X'100)
						// X'100
						qu2.cursor = new_que_e2_n(que);
						qu2.key = qu.key + query[q].bd[query[q].idx[2]] - query[q].bd[query[q].idx[1]];
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ MASK(query[q].idx[2])) ^ MASK(query[q].idx[1]);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
						// X11
						qu.key = qu.key + query[q].bd[query[q].idx[0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ MASK(query[q].idx[0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						#ifdef USE_TOP
						enq_c2_n_after_deq(&qu2, que);
						#else
						enq_c2_n(&qu2, que);
						#endif
						enq_c2_n(&qu, que);
					} else {
						qu.key = qu.key + query[q].bd[query[q].idx[0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ MASK(query[q].idx[0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						#ifdef USE_TOP
						enq_c2_n_after_deq(&qu, que);
						#else
						enq_c2_n(&qu, que);
						#endif
					}
					break;
				case 3: // X11 -> nothing to do
					#ifdef CURSOR_STACK
					que->cursor_stack[++que->stack_top] = qu.cursor;
//					printf("returned cursor = %d\n", qu.cursor);
					#endif
					#ifdef USE_TOP
					deq_c2_n_del(que);
					#endif
					break;
				}
				#endif
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
//				num_enumerated += enumerated;
//				num_remained += que->qsize;
			}
END_SEARCH: ;
		}
	}
//	printf("max detail size = %d\n", max_detail_size); // exit(0);
//	printf("enumerated = %d, remained = %d, shurinked = %6.2lf%%\n", num_enumerated, num_remained, (double)(num_enumerated - num_remained) / num_enumerated); exit(0);
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// カーソルを使った優先度付き待ち行列を使用するimamuraバージョン（delayed enq を用いている）
// 一部は，樋口のideaによる
	// Heap をうまく利用して，stackを使わないで，再利用を行う (shino)
double precision_resample_cursor_2_s(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif
//	int num_enumerated = 0, num_remained = 0;

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	static struct_que_c2_s *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que_c2_s *)malloc(sizeof(struct_que_c2_s) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
		for(q = 0; q < 16; q++) {
			make_init_que_c2_s(&que_pool[q], QSIZE + 1);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		register sketch_type s;
//		QUE_c2 *qu, *qu2;
		dist_type key;
		sketch_type pt;
		struct_que_c2_s *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均
//		int enumerated = 0;

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
			
			s = query[q].sketch;

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j; 
					break;
				}
			}
			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
				if(k >= trunc_K) goto END_SEARCH;
			}

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
			{
				for(tn = 0; tn < 16; tn++) {
					if(used_q[tn] == 0) {
						que = &que_pool[tn];
						used_q[tn] = 1;
						break;
					}
					if(tn == 16) {
						fprintf(stderr, "No free que in que_pool\n");
						exit(0);
					}
				}
			}

			make_empty_que_c2_s(que);

			// enq pattern of 0...10
//			qu = new_que_e2_s(que);
//			qu->key = query[q].bd[query[q].idx[1]];
//			que->details[qu->cursor].sk = query[q].sketch ^ (1 << query[q].idx[1]);
//			que->details[qu->cursor].pt = 1 << 1; // pt = "0...00000010"
			que->element[0].key = query[q].bd[query[q].idx[1]];
			que->details[que->element[0].cursor].sk = query[q].sketch ^ (1 << query[q].idx[1]);
			que->details[que->element[0].cursor].pt = 1 << 1; // pt = "0...00000010"
			enq_c2_s(que);

//			while(deq_c2_s(&qu, que)) {
			while(que->qsize) {
//	START_SEARCH:
//				enumerated++;
//				printf("cursor = %d, pt = %d, key = %d, qsize = %d\n", qu.cursor, que->details[qu.cursor].pt, qu.key, que->qsize);
				key = que->element[0].key;
				if(key > score) {
					if(k >= trunc_K) goto FREE_Q;
					score = key;
				}
				s = que->details[que->element[0].cursor].sk;
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}
				pt = que->details[que->element[0].cursor].pt;
				switch(pt & 3) {
				case 0: // X00 -> enq(X01) and enq(Y10^{m+1}) if X00 = Y010^m
					{
						int m = lsb_pos(pt);
						if(m > 0 && m < PJT_DIM - 1 && !(pt & MASK(m + 1))) {
							// Y010^m -> Y10^{m+1}
							que->element[0].key = key + query[q].bd[query[q].idx[m + 1]] - query[q].bd[query[q].idx[m]];
							que->details[que->element[0].cursor].sk = (s ^ MASK(query[q].idx[m + 1])) ^ MASK(query[q].idx[m]);
							que->details[que->element[0].cursor].pt = pt + MASK(m);
							// Y010^m -> Y010^{m-1}1
							que->element[que->qsize].key = key + query[q].bd[query[q].idx[0]];
							que->details[que->element[que->qsize].cursor].sk = s ^ MASK(query[q].idx[0]);
							que->details[que->element[que->qsize].cursor].pt = pt + 1;
							enq_c2_s_after_deq(que);
							enq_c2_s(que);
						} else {
							que->element[0].key = key + query[q].bd[query[q].idx[0]];
							que->details[que->element[0].cursor].sk = s ^ MASK(query[q].idx[0]);
							que->details[que->element[0].cursor].pt = pt + 1;
							enq_c2_s_after_deq(que);
						}
					}
					break;
				case 1: // X01 -> enq(X10) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
					que->element[0].key = key + query[q].bd[query[q].idx[1]] - query[q].bd[query[q].idx[0]];
					que->details[que->element[0].cursor].sk = (s ^ MASK(query[q].idx[1])) ^ MASK(query[q].idx[0]);
					que->details[que->element[0].cursor].pt = pt + 1;
					enq_c2_s_after_deq(que);
					break;
				case 2: // X10 -> enq(X11) and (enq(x100) only if X10 = X'010)
					if(!(pt & 4)) { // X'010 -> enq(X'100)
						// X'100
						que->element[0].key = key + query[q].bd[query[q].idx[2]] - query[q].bd[query[q].idx[1]];
						que->details[que->element[0].cursor].sk = (s ^ MASK(query[q].idx[2])) ^ MASK(query[q].idx[1]);
						que->details[que->element[0].cursor].pt = pt + 2;
						// X11
						que->element[que->qsize].key = key + query[q].bd[query[q].idx[0]];
						que->details[que->element[que->qsize].cursor].sk = s ^ MASK(query[q].idx[0]);
						que->details[que->element[que->qsize].cursor].pt = pt + 1;
						enq_c2_s_after_deq(que);
						enq_c2_s(que);
					} else {
						que->element[0].key = key + query[q].bd[query[q].idx[0]];
						que->details[que->element[0].cursor].sk = s ^ MASK(query[q].idx[0]);
						que->details[que->element[0].cursor].pt = pt + 1;
						enq_c2_s_after_deq(que);
					}
					break;
				case 3: // X11 -> nothing to do
					deq_c2_s_del(que);
					break;
				}
			}
FREE_Q:
//			make_init_que_c2_s(que, QSIZE + 1 /* 50000 enumerated + 10 */);
			#pragma omp critical 
			{
				used_q[tn] = 0;
//				num_enumerated += enumerated;
//				num_remained += que->qsize;
			}
END_SEARCH: ;
		}
	}
//	printf("max detail size = %d\n", max_detail_size); // exit(0);
//	printf("enumerated = %d, remained = %d, shurinked = %6.2lf%%\n", num_enumerated, num_remained, (double)(num_enumerated - num_remained) / num_enumerated); exit(0);
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// 一部は，樋口のideaによる
// Heap を stack として使って，再利用を行う (shino)
double precision_resample_s(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif
//	int num_enumerated = 0, num_remained = 0;

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	static struct_que_s *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que_s *)malloc(sizeof(struct_que_s) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
		for(q = 0; q < 16; q++) {
			make_empty_que_s(&que_pool[q]);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		register sketch_type s;
//		QUE_c2 *qu, *qu2;
		dist_type key;
		sketch_type pt;
		struct_que_s *que = NULL;
		int *element;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均
//		int enumerated = 0;

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
			
			s = query[q].sketch;

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j; 
					break;
				}
			}
			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
				if(k >= trunc_K) goto END_SEARCH;
			}

			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
			{
				for(tn = 0; tn < 16; tn++) {
					if(used_q[tn] == 0) {
						que = &que_pool[tn];
						used_q[tn] = 1;
						break;
					}
					if(tn == 16) {
						fprintf(stderr, "No free que in que_pool\n");
						exit(0);
					}
				}
			}

			make_init_que_s(que, 100000);
			element = que->element;

			// enq pattern of 0...10
			KEY(element[0]) = query[q].bd[query[q].idx[1]];
			SK(element[0]) = query[q].sketch ^ (1 << query[q].idx[1]);
			PT(element[0]) = 1 << 1; // pt = "0...00000010"
			enq_s(que);

//			while(deq_c2_s(&qu, que)) {
			while(que->qsize) {
//	START_SEARCH:
//				enumerated++;
//				printf("cursor = %d, pt = %d, key = %d, qsize = %d\n", qu.cursor, que->details[qu.cursor].pt, qu.key, que->qsize);
				key = KEY(element[0]);
				if(key > score) {
					if(k >= trunc_K) goto FREE_Q;
					score = key;
				}
				s = SK(element[0]);
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[s]; j < bkt[s + 1]; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}
				pt = PT(element[0]);
				switch(pt & 3) {	// X01 or X11
					int m;
				case 3:	// X11
					deq_s_del(que);
					break;
				case 1:	// X01 -> enq_after_deq(X10)
					KEY(element[0]) = key + query[q].bd[query[q].idx[1]] - query[q].bd[query[q].idx[0]];
					SK(element[0]) = (s ^ MASK(query[q].idx[1])) ^ MASK(query[q].idx[0]);
					PT(element[0]) = pt + 1;
					enq_s_after_deq(que);
					break;
				case 0: case 2:	// X00 or X10
					m = lsb_pos(pt);
					if(m < PJT_DIM - 1 && !(pt & MASK(m + 1))) {
						// Y010^m -> Y10^{m+1}
						KEY(element[que->qsize]) = key + query[q].bd[query[q].idx[m + 1]] - query[q].bd[query[q].idx[m]];
						SK(element[que->qsize]) = (s ^ MASK(query[q].idx[m + 1])) ^ MASK(query[q].idx[m]);
						PT(element[que->qsize]) = pt + MASK(m);
						enq_s(que);
					}
					KEY(element[0]) = key + query[q].bd[query[q].idx[0]];
					SK(element[0]) = s ^ MASK(query[q].idx[0]);
					PT(element[0]) = pt + 1;
					enq_s_after_deq(que);
					break;
				}
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
//				num_enumerated += enumerated;
//				num_remained += que->qsize;
			}
END_SEARCH: ;
		}
	}
//	printf("max detail size = %d\n", max_detail_size); // exit(0);
//	printf("enumerated = %d, remained = %d, shurinked = %6.2lf%%\n", num_enumerated, num_remained, (double)(num_enumerated - num_remained) / num_enumerated); exit(0);
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// スケッチ（あるいは質問スケッチとのXOR)を構造体データに記録しないで，格納カーソルで表現する
double precision_resample_cursor_sk(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
#if SKETCH_PARTITION == 3
	int pruned = 0, total_pruned = 0;
#endif

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;

	static struct_que_sk *que_pool = NULL;
	int used_q[16] = {0};
	if(que_pool == NULL) {
		que_pool = (struct_que_sk *)malloc(sizeof(struct_que_c) * 16);
		if(que_pool == NULL) {
			fprintf(stderr, "Malloc for que failed.\n");
			exit(0);
		}
		for(q = 0; q < 16; q++) {
			make_empty_que_sk(&que_pool[q]);
		}
	}

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		sketch_type qu, qu2, qu3;
		struct_que_sk *que = NULL;
		int k;
		dist_type score = 0;
		int mid_k; // 同じスケッチを持つデータが現れたときの初めの k と最後の k の平均

		if(lb == 1 || lb >= 5) { // L1
			k = 0;
#if SKETCH_PARTITION == 3
			pruned = 0;
#endif
			s = query[q].sketch;
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1] && k < trunc_K; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; // k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j; 
					break;
				}
			}
			s = s ^ (1 << query[q].idx[0]);
			if(query[q].bd[query[q].idx[0]] > score) {
				score = query[q].bd[query[q].idx[0]];
			}
			mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
			for(int j = bkt[s]; j < bkt[s + 1] && k < trunc_K; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
					goto END_SEARCH; //k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}

			int tn;
			#pragma omp critical 
			{
				for(tn = 0; tn < 16; tn++) {
					if(used_q[tn] == 0) {
						que = &que_pool[tn];
						used_q[tn] = 1;
						break;
					}
					if(tn == 16) {
						fprintf(stderr, "No free que in que_pool\n");
						exit(0);
					}
				}
			}

			make_empty_que_sk(que);

			qu = 1 << query[q].idx[1];
			que->key[qu] = query[q].bd[query[q].idx[1]];
			que->idx[qu] = 2;
			enq_sk(qu, que);
			while(deq_sk(&qu, que) && k < trunc_K) {
				s = query[q].sketch ^ qu;
				if(que->key[qu] > score) {
					score = que->key[qu];
				}
				mid_k = k + (bkt[s + 1] - bkt[s]) / 2; // 最初の k と sk_num[s] 後の k + sk_num[s] の平均
				for(int j = bkt[s]; j < bkt[s + 1] && k < trunc_K; j++, k++) {
					if(idx[j] == query[q].nearest) {
						num_success++;
						sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
						ans[q].idx = q + 1;
						ans[q].dist = k + 1;
						goto FREE_Q; //k = db_header->data_num; // 検索終了
					} else if(idx[j] > query[q].nearest) {
						k += bkt[s + 1] - j;
						break;
					}
				}
				if(que->idx[qu] < PJT_DIM) {
					qu2 = (qu ^ (1 << query[q].idx[que->idx[qu] - 1])) ^ (1 << query[q].idx[que->idx[qu]]);
					que->key[qu2] = que->key[qu] + query[q].bd[query[q].idx[que->idx[qu]]] - query[q].bd[query[q].idx[que->idx[qu] - 1]];
					que->idx[qu2] = que->idx[qu] + 1;
					enq_sk(qu2, que);
					if(!(qu & (1 << query[q].idx[que->idx[qu] - 2]))) {
						qu3 = qu ^ (1 << query[q].idx[que->idx[qu] - 2]);
						que->key[qu3] = que->key[qu] + query[q].bd[query[q].idx[que->idx[qu] - 2]];
						que->idx[qu3] = que->idx[qu];
						enq_sk(qu3, que);
					}
				}
			}
FREE_Q:
			#pragma omp critical 
			{
				used_q[tn] = 0;
			}
END_SEARCH: ;
		}
	}
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

// 再サンプリング(大きさ = resampling_size)の質問に対する精度（score_infを用いる）
double precision_resample_inf(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合
	//                         1 → 重み付きスコア合計
	int i, q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;

	prepare_bucket_db(2, sketch, db_header->data_num, &idx, &bkt);
//	fprintf(stderr, "trunc_K = %d, average number of elements in nonempty buckets = %d\n", trunc_K, nonempty_ave);

	#pragma omp parallel for
	for(i = 0; i < resampling_size; i++)
		ans[i].dist = UINT_MAX;


	#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: num_success, sum) schedule(static)
	#endif
	for(q = 0; q < resampling_size; q++) {
		// local variables for parallel
		sketch_type s;
		bd_type bd_idx[PJT_DIM];
		int i, j, k;

		bound_with_idx(query[q].data, bd_idx);
		qsort(bd_idx, PJT_DIM, sizeof(bd_type), comp_bd_type);

		s = query[q].sketch;
//		printf("(a) s = "); print_bin(s); printf(", sk_num = %d\n", sk_num[s]);
		k = 0;
		for(j = bkt[s]; j < bkt[s + 1] && k < trunc_K; j++, k++) {
			if(idx[j] == query[q].nearest) {
				num_success++;
				sum += log(SCORE_FACTOR * db_header->data_num / (k + 1));
				ans[q].idx = q + 1;
				ans[q].dist = k + 1;
//				fprintf(stderr, "(a) found k = %d\n", k);
				k = db_header->data_num; // 検索終了
			} else if(idx[j] > query[q].nearest) {
				k += bkt[s + 1] - j;
				break;
			}
		}
//		fprintf(stderr, "(a) OK\n");
		
		for(i = 0; i < BIT && k < trunc_K; i++) {
			s ^= 1 << bd_idx[bitcnt_tbl[i ^ (i + 1)] - 1].idx;	
//			printf("(b) s = "); print_bin(s); printf(", sk_num = %d, k = %6d\n", sk_num[s], k);
			for(j = bkt[s]; j < bkt[s + 1] && k < trunc_K; j++, k++) {
				if(idx[j] == query[q].nearest) {
					num_success++;
					sum += log(SCORE_FACTOR * db_header->data_num / (k + 1));
					ans[q].idx = q + 1;
					ans[q].dist = k + 1;
//					fprintf(stderr, "(b) found k = %d\n", k);
					k = db_header->data_num; // 検索終了
				} else if(idx[j] > query[q].nearest) {
					k += bkt[s + 1] - j;
					break;
				}
			}
		}
	}
	*prec = (double)100 * num_success / resampling_size;
	return (mode == 1 ? sum / resampling_size: (double)100 * num_success / resampling_size);
}

#ifdef EVAL_BY_SEQUENTIAL_FILTERING
// AIRやLSによるピボット探索時に用いる Sequential Filtering に基づく，精度（recall）の評価		
double precision_resample_by_sequential_filtering(int obj, int resampling_size, int mode, double *prec) {
	// mode: 評価の仕方の指定：0 → 候補数割合を K(%) = obj / 10 としたときの正解数の割合（未実装）
	//                         1 → 重み付きスコア合計
	int q;
	int trunc_K = (double)obj / 1000 * db_header->data_num;
	int num_success = 0;
	double sum = 0;
// #if SKETCH_PARTITION == 3 // 当面無視
//	int pruned = 0, total_pruned = 0;
// #endif

	#pragma omp parallel for
	for(q = 0; q < resampling_size; q++)
		ans[q].dist = UINT_MAX;

#ifdef TRAINING_PARALLEL
	#pragma omp parallel for reduction (+: sum) schedule(static)
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
		for( ; i < db_header->data_num / 8; i++) {
			score = priority(sketch[i], &query[q]);
			if(score < query[q].answer_score) {
				n_low++;
			} else if(score == query[q].answer_score) {
				n_eq++;
			}
		}
		#ifdef TRUNC_K
		if(n_low > trunc_K) continue;
		#endif
		for( ; i < db_header->data_num; i++) {
			score = priority(sketch[i], &query[q]);
			if(score < query[q].answer_score) {
				n_low++;
			} else if(score == query[q].answer_score) {
				n_eq++;
			}
		}
		#ifdef TRUNC_K
		if(n_low > trunc_K) continue;
		#endif
		mid_k = n_low + n_eq / 2;
		if(mid_k < trunc_K) {
			num_success++;
		}
		#ifdef TRUNC_K
		sum += log(SCORE_FACTOR * db_header->data_num / (mid_k + 1));
		#else
		sum += log(db_header->data_num / (mid_k + 1));
		#endif
		ans[q].idx = q + 1;
		ans[q].dist = mid_k + 1;
	}

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
