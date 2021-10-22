#!/bin/bash
if [ $# != 5 ] && [ $# != 4 ] ; then
echo "Usage> スクリプト.sh TRIAL_LS WIDTH SEED RETRY K"
echo "1. TRIAL_AIR  = AIRにおける各ピボット選択の試行回数 (TRIAL2)"
echo "2. WIDTH      = スケッチの幅（ビット数）" 
echo "3. SEED       = 初期設定のSEED"
echo "4. RETRY      = QBP による候補選択数（REPLACE_TRIAL）"
echo "5. K (optional) = 評価時の打切り候補数の割合 ％ の10倍"
exit 1
fi

if [ $# == 5 ] ; then
tk=$5
else
tk=10
fi

db_dir="/mnt/f/DISA_h/dataset"
qr_dir="/mnt/g/DISA_h/query"
pv_dir="../pivot"
pr_dir=/mnt/g/DISA_h/program/sketch

sd="$db_dir/fc6_00_one_32th.ftr"
q1="$qr_dir/fc6_q_00_04.ftr"
a1="$qr_dir/fc6_q_00_04_one_32th.csv"

if [ $# == 5 ] ; then
sk="SK/pivot_AIR"$1"_w"$2"_s"$3"_one_32th_k"${tk}"_seq.sk"
else
sk="SK/pivot_AIR"$1"_w"$2"_s"$3"_one_32th_no_trunc_K_seq.sk"
fi
#pv="SK/pivot_AIR"$1"_w"$2"_s"$3"_16th.pivot.csv"
r1="result_AIR_w${2}.csv"

cflags="-O3 -fopenmp -Wall -Wno-strict-overflow"
#cflags="-O3 -Wall -Wno-strict-overflow"
#cflags="$cflags -DEVAL_PRECISION=precision_resample_inf"
#cflags="$cflags -DEVAL_PRECISION=precision_resample_cursor_2_n"
cflags="$cflags -DEVAL_PRECISION=precision_resample_by_sequential_filtering"
cflags="$cflags -DEVAL_BY_SEQUENTIAL_FILTERING"
cflags="$cflags -DENQ_P=enq_p2 -DDEQ_P=deq_p2"
cflags="$cflags -DREPLACE_WHOLE -DREPLACE_TRIAL=$4"
cflags="$cflags -DTRAINING_PARALLEL -DTEST_PARALLEL"
#cflags="$cflags -DMEASURE_SEARCH_TIME"
cflags="$cflags -DWITHOUT_MATCHING -DUSE_SORTED_DB=0 -DVAR_FLIP -DT_POWER=1"
cflags="$cflags -DFTR_DIM=4096 -DPJT_DIM=${2}"
cflags="$cflags -DEVAL_MODE=1"
cflags="$cflags -DSKETCH_PARTITION=2"
cflags="$cflags -DNUM_TRIAL1=10"
cflags="$cflags -DN0=10 -DNUM_TRIAL2=$1"
# NUM_TRIAL2 は AIR の試行回数（0 にすると LS のみになる）
cflags="$cflags -DREUSE=10"
cflags="$cflags -DNUM_LS=0"
# NUM_LS は LS の試行回数（0 にすると AIR になる．両方指定すると，AIR ＋ LS になる）
cflags="$cflags -DFLIP_NUM=10"
cflags="$cflags -DSCORE_FACTOR=0.4"
if [ $# == 5 ] ; then
cflags="$cflags -DTRUNC_K"
fi
echo $cflags

grep "model name" /proc/cpuinfo | tail -1

cp -v -u -p $pr_dir/*.c $pr_dir/*.h ./program
gcc $cflags ./program/make_Sketch_by_AIR.c -o make_Sketch_by_AIR -lm

if [ $? == 1 ] ; then
exit 1
fi

grep "model name" /proc/cpuinfo | tail -1
echo "replace trial = " $4
#                                              score_1 partition 10*K seed
time ./make_Sketch_by_AIR $sd $q1 $a1 $sk       1         2      $tk   $3

set +x
exit 0

