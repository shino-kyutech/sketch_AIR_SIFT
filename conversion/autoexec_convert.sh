#!/bin/bash

gcc -O3 -DBUFFER_SIZE=100000 -fopenmp convert_Deep1B_to_FTR_DB_and_Q.c -o convert_Deep1B_to_FTR_DB_and_Q

fbin=/mnt/poseidon/share2/Deep1B/base.1B.fbin
for n in {0..99} ; do
if [ $n -lt 10 ] ; then
ftr=/mnt/e/Deep1B/base_0$n.ftr
else
ftr=/mnt/e/Deep1B/base_$n.ftr
fi
num=10000000
from=$(($n * ${num}))
#./convert_Deep1B_to_FTR_DB_and_Q /mnt/poseidon/share2/Deep1B/base.10M.fbin base_10M.ftr 0 0
./convert_Deep1B_to_FTR_DB_and_Q ${fbin} ${ftr} ${from} ${num}
#echo ftr = ${ftr}, from = ${from}, num = ${num}
done

