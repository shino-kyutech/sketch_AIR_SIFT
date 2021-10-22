#!/bin/bash
# 入力:		Deep1B のメインの10Mデータセット"base.10M.fbin"（1000万個）
# 出力:		10個の部分データセット（100万個），base10M_00.ftr, ... , base10M_09.ftr
# Data_id：	1350000000 から 1359999999（オフセット無し）

gcc -O3 -DBUFFER_SIZE=100000 -fopenmp convert_Deep1B_to_FTR_DB_and_Q.c -o convert_Deep1B_to_FTR_DB_and_Q
if [ $? == 1 ] ; then
echo compile error
exit 1
fi

fbin=/mnt/poseidon/share2/Deep1B/base.10M.fbin
for n in {0..9} ; do
if [ $n -lt 10 ] ; then
ftr=/mnt/e/Deep1B/dataset10M/base10M_0$n.ftr
else
ftr=/mnt/e/Deep1B/dataset10M/base10M_$n.ftr
fi
num=1000000
from=$(($n * ${num}))
offset=1350000000
./convert_Deep1B_to_FTR_DB_and_Q ${fbin} ${ftr} ${from} ${num} ${offset}
done

