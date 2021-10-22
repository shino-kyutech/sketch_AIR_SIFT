#!/bin/bash
# 入力:		Deep1B の学習用データセット"learn.350M.fbin"（3.5億個）
# 出力:		35個の学習用部分データセット（1000万個），learn_00.ftr, ... , learn_34.ftr
# Data_id：	1000000000 から 1349999999（オフセット = 1000000000）

gcc -O3 -DBUFFER_SIZE=100000 -fopenmp convert_Deep1B_to_FTR_DB_and_Q.c -o convert_Deep1B_to_FTR_DB_and_Q
if [ $? == 1 ] ; then
echo compile error
exit 1
fi

fbin=/mnt/poseidon/share2/Deep1B/learn.350M.fbin
for n in {0..34} ; do
if [ $n -lt 10 ] ; then
ftr=/mnt/e/Deep1B/learn/learn_0$n.ftr
else
ftr=/mnt/e/Deep1B/learn/learn_$n.ftr
fi
num=10000000
from=$(($n * ${num}))
offset=1000000000
./convert_Deep1B_to_FTR_DB_and_Q ${fbin} ${ftr} ${from} ${num} ${offset}
done
