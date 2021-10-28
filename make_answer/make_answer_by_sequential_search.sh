#!/bin/bash
if [ $# -lt 5 ] ; then
echo "Usage> スクリプト.sh query.ftr num_k #threads storage_type mode dataset_1.ftr ... dataset_n.ftr"
echo "1. query.ftr   = query file" 
echo "2. num_k       = k"
echo "3. #threads    = number of threads"
echo "4. mode        = PRINTOUT (to print search results) or NO_PRINTOUT"
echo "5 - dataset_1.ftr, ... , dataset_n file = dataset files" 
exit 1
fi

ds_dir="/mnt/e/Deep1B/dataset"
qr_dir="/mnt/e/Deep1B/query"
pr_dir="../src"

qr="$qr_dir/$1"; shift
if [ ! -e $qr ]; then
  echo query file = $qr does not exist.
  exit
fi

nk="$1"; shift
nt="$1"; shift

if [ $nt == 1 ] ; then
cflags="-O3 -Wall -Wno-strict-overflow"
else
cflags="-O3 -fopenmp -Wall -Wno-strict-overflow"
cflags="$cflags -DNUM_THREADS=$nt"
fi

cflags="$cflags -DCOMPILE_TIME"
cflags="$cflags -DDEEP1B"
cflags="$cflags -DFTR_ON_SECONDARY_MEMORY"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"
cflags="$cflags -DWITHOUT_FILTERING"
cflags="$cflags -DPARTITION_TYPE_QBP"
cflags="$cflags -DANSWER_WITH_DATA_ID"
cflags="$cflags -DQUERY=\"$qr\""
cflags="$cflags -DNUM_K=$nk"

if [ $1 == "PRINTOUT" ] ; then
cflags="$cflags -DPRINTOUT"
fi
shift

echo $cflags

ds=" "
count=1
while [ "$#" -ge "1" ]; do
  if [ ! -e $ds_dir/$1 ]; then
    echo dataset file = $ds_dir/$1 does not exist
    exit 1
  fi
  ds="$ds $ds_dir/$1 "
  shift
  let count=$count+1
done

echo $ds

gcc $cflags $pr_dir/make_answer_by_sequential_search.c $pr_dir/ftr.c $pr_dir/kNN_search.c -o make_answer_by_sequential_search

if [ $? == 1 ] ; then
echo "compile error!"
exit 1
fi

time ./make_answer_by_sequential_search $ds

