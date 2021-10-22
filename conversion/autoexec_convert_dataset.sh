#!/bin/bash
# 入力:		Deep1B のメインのデータセット"base.1B.fbin"（10億個）
# 出力:		100個の部分データセット（1000万個），base_00.ftr, ... , base_99.ftr
# Data_id：	0 から 999999999（オフセット無し）

gcc -O3 -DBUFFER_SIZE=100000 -fopenmp convert_Deep1B_to_FTR_DB_and_Q.c -o convert_Deep1B_to_FTR_DB_and_Q
if [ $? == 1 ] ; then
echo compile error
exit 1
fi

fbin=/mnt/poseidon/share2/Deep1B/base.1B.fbin
for n in {0..99} ; do
if [ $n -lt 10 ] ; then
ftr=/mnt/e/Deep1B/dataset/base_0$n.ftr
else
ftr=/mnt/e/Deep1B/dataset/base_$n.ftr
fi
num=10000000
from=$(($n * ${num}))
offset=0
./convert_Deep1B_to_FTR_DB_and_Q ${fbin} ${ftr} ${from} ${num} ${offset}
done

