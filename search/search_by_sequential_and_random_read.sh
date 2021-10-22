#!/bin/bash
if [ $# -ne 12 ] ; then
echo "Usage>"
echo "1.  query      = query file (without .ftr)" 
echo "2.  from       = file number of FROM"
echo "3.  to         = file number of TO"
echo "4.  ftr_on     = RAM -> main memory, SSD -> secondary memory (SSD), HDD -> secondary memory (HDD)"
echo "5.  num_d      = number of data to be searched"
echo "6.  num_k      = k (number of neighbors)"
echo "7.  num_q      = number of queries (0 -> all)"
echo "8.  repeat     = number of repetition"
echo "9.  block_size = number of records to read together" 
echo "10. #threads   = number of threads for parallel process"
echo "11. mode       = SEQUENTIAL -> sequelntial read, RANDOM -> random read"
echo "12. seed       = random seed"
exit 1
fi

ulimit -Sn 4000

#set -x

query=$1; shift
from=$1; shift
to=$1;shift
#range=$1; shift
fm=$1; shift
numd=$1; shift
numk=$1; shift
numq=$1; shift
repeat=$1; shift
blks=$1; shift
numt=$1; shift
mode=$1; shift
seed=$1; shift

if [ $fm == RAM ] ; then
source ../src/set_directory.sh SSD
ftron="MAIN_MEMORY"
elif [ $fm == SSD ] ; then
source ../src/set_directory.sh INTEL
ftron="SECONDARY_MEMORY"
elif [ $fm == HDD ] ; then
source ../src/set_directory.sh HDD
ftron="SECONDARY_MEMORY"
else
echo "invalid parameter: ftr_on = RAM -> main memory, SSD -> secondary memory (SSD), HDD -> secondary memory (HDD)"
exit
fi

if [ $mode != SEQUENTIAL ] && [ $mode != RANDOM ] && [ $mode != BOTH ] ; then
echo "invalid access mode: SEQUENTIAL -> sequential read, RANDOM -> random read, BOTH -> sequential and random"
exit
fi

if [ $from == ALL ] ; then
ds="${ds_dir}/base_00_99.ftr"
elif [ $from == w20 ]; then
ds="${sf_dir}/pivot_QBP_t1000_w20_sd00_ss10000_seed1_00_99.sftr"
elif [ $from == w22 ]; then
ds="${sf_dir}/pivot_QBP_t1000_w22_sd00_ss10000_seed1_00_99.sftr"
elif [ $from == w24 ] ; then
ds="${sf_dir}/pivot_QBP_t1000_w24_sd00_ss10000_00_99.sftr"
elif [ $from == w26 ] ; then
ds="${sf_dir}/pivot_QBP_t1000_w26_sd00_ss10000_00_99.sftr"
else
range=$(../src/expand_fn.sh "${ds_dir}/base_" ${from} ${to} ".ftr" 0)
ds=$(../src/expand_fn.sh "${ds_dir}/base_" ${from} ${to} ".ftr" 1)
fi

for f in $ds ; do
if [ ! -e $f ] ; then
  echo ftr file = $f does not exist.
  exit
fi
done

qr="$qr_dir/${query}.ftr"
if [ ! -e $qr ]; then
  echo query file = $qr does not exist.
  exit
fi

if [ $numt == 1 ] ; then
cflags="-O3 -Wall -Wno-strict-overflow"
else
cflags="-O3 -fopenmp -Wall -Wno-strict-overflow"
cflags="$cflags -DNUM_THREADS=$numt"
fi

cflags="$cflags -DCOMPILE_TIME"
cflags="$cflags -DDEEP1B"
cflags="$cflags -DWITHOUT_FILTERING"
#cflags="$cflags -DFTR_DIM=4096"
#cflags="$cflags -DPJT_DIM=16"
cflags="$cflags -DPARTITION_TYPE_QBP"
cflags="$cflags -DFTR_ON_${ftron}"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"
cflags="$cflags -DBLOCK_SIZE=${blks}"
cflags="$cflags -DNUM_D=$numd"
cflags="$cflags -DNUM_K=$numk"
cflags="$cflags -DNUM_Q=$numq"
cflags="$cflags -DREPEAT=$repeat"
cflags="$cflags -DQUERY_FILE=\"$qr\""
cflags="$cflags -DACCESS_MODE=\"${mode}\""
cflags="$cflags -DSEED=$seed"

echo $cflags

l1=bit_op
l2=ftr
l3=kNN_search
l4=e_time
#l5=quick
#l6=bucket

lp="$pr_dir/$l1.c $pr_dir/$l2.c $pr_dir/$l3.c $pr_dir/$l4.c"
#lp="$pr_dir/$l1.c $pr_dir/$l2.c $pr_dir/$l3.c $pr_dir/$l4.c $pr_dir/$l5.c"

gcc $cflags -c $lp

if [ $? == 1 ] ; then
echo compile error for libraly
exit 1
fi

lb="$l1.o $l2.o $l3.o $l4.o -lm"
#lb="$pr_dir/$l1.o $pr_dir/$l2.o $pr_dir/$l3.o $pr_dir/$l4.o $pr_dir/$l5.o"

pr=search_by_sequential_and_random_read

gcc $cflags $pr_dir/$pr.c $lb -o $pr

if [ $? == 1 ] ; then
exit 1
fi

time ./$pr $ds

