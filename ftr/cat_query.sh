#!/bin/bash

out=$1; shift
from=$1; shift
to=$1; shift

cflags="-O3 -Wall -Wno-strict-overflow"
cflags="$cflags -DCOMPILE_TIME"
cflags="$cflags -DDEEP1B"
cflags="$cflags -DFTR_ON_SECONDARY_MEMORY"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"
cflags="$cflags -DSEQUENTIAL_FILTERING"
cflags="$cflags -DPARTITION_TYPE_GHP"


gcc $cflags ../src/cat_FTR.c ../src/ftr.c -o cat_FTR
if [ $? == 1 ] ; then
exit 1
fi

ds=$(../src/expand_fn.sh "/mnt/e/Deep1B/query/query_" ${from} ${to} ".ftr" 1)

./cat_FTR $out $ds

