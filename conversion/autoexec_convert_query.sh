#!/bin/bash
# 入力:		Deep1B の質問用データセット"query.public.10K.fbin"（1万個）
# 出力:		10個の学習用部分データセット（1000個），query_00.ftr, ... , query_09.ftr
# Data_id：	1500000000 から 1500009999（オフセット = 1500000000）

gcc -O3 -DBUFFER_SIZE=1000 -fopenmp convert_Deep1B_to_FTR_DB_and_Q.c -o convert_Deep1B_to_FTR_DB_and_Q
if [ $? == 1 ] ; then
echo compile error
exit 1
fi

fbin=/mnt/poseidon/share2/Deep1B/query.public.10K.fbin
for n in {0..9} ; do
if [ $n -lt 10 ] ; then
ftr=/mnt/e/Deep1B/query/query_0$n.ftr
else
ftr=/mnt/e/Deep1B/query/query_0$n.ftr
fi
num=1000
from=$(($n * ${num}))
offset=1500000000
./convert_Deep1B_to_FTR_DB_and_Q ${fbin} ${ftr} ${from} ${num} ${offset}
done

