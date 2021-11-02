#!/bin/bash
if [ $# -ne 3 ] ; then
echo "Usage>"
echo "1. pivot       = pivot file"
echo "2. range       = range of ftr files" 
echo "3. block_size  = number of data in block access"
exit 1
fi

#set -x

ulimit -Sn 4000

source ../src/set_directory.sh HDD

pivot=$1; shift

pv="$pv_dir/${pivot}.csv"
if [ ! -e $pv ]; then
  echo pivot file = $pv does not exist.
  exit
fi
wd=$(../src/pivot_property.sh -w $pv)
dm=$(../src/pivot_property.sh -d $pv)
pt=$(../src/pivot_property.sh -p $pv)
np=$(../src/pivot_property.sh -n $pv)

range=$1; shift
bk="$bk_dir/${pivot}_${range}.bkt"
if [ ! -e $bk ]; then
  echo bucket file = $bk does not exist.
  exit
fi

block=$1; shift

sf="$sf_dir/${pivot}_${range}.sftr"

#cflags0="-O3 -Wall -Wno-strict-overflow"
cflags0="-O3 -fopenmp -Wall -Wno-strict-overflow"
cflags0="$cflags0 -DNUM_THREADS=16"

cflags0="$cflags0 -DCOMPILE_TIME"
cflags0="$cflags0 -DDEEP1B"
cflags0="$cflags0 -DFTR_DIM=$dm"
cflags0="$cflags0 -DPJT_DIM=$wd"

if [ $wd -gt 64 ] ; then
cflags0="$cflags0 -DEXPANDED_SKETCH"
elif [ $wd -gt 32 ] ; then
cflags0="$cflags0 -DWIDE_SKETCH"
else
cflags0="$cflags0 -DNARROW_SKETCH"
fi

if [ $pt == 3 ] ; then
cflags0="$cflags0 -DPARTITION_TYPE_PQBP"
cflags0="$cflags0 -DNUM_PART=$np"
else
cflags0="$cflags0 -DPARTITION_TYPE_QBP"
fi

cflags0="$cflags0 -DBLOCK_SIZE=${block}"
cflags="$cflags -DFTR_ON_SECONDARY_MEMORY"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"
#cflags="$cflags -DBLOCK_SIZE=1"
cflags="$cflags -DSEQUENTIAL_FILTERING"
cflags="$cflags -DBUCKET_FILE=\"$bk\""
cflags="$cflags -DSFTR_FILE=\"$sf\""

echo $cflags

ds=$(../src/expand_filenames.sh "${ds_dir}/base_" $range ".ftr")

echo $ds

l1=bit_op
l2=ftr
l3=kNN_search
l4=sketch
l5=quick

lp="$pr_dir/$l1.c $pr_dir/$l2.c $pr_dir/$l3.c $pr_dir/$l4.c $pr_dir/$l5.c"

gcc $cflags0 $cflags -c $lp

if [ $? == 1 ] ; then
echo compile error for libraly
exit 1
fi

lb="$l1.o $l2.o $l3.o $l4.o $l5.o -lm"

pr=make_sorted_ftr

gcc $cflags0 $cflags $pr_dir/$pr.c $lb -o $pr

if [ $? == 1 ] ; then
exit 1
fi

time ./$pr $ds

