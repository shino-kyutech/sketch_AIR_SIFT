#!/bin/bash
if [ $# -ne 6 ] ; then
echo "Usage>"
echo "1. ratio      = ratio of shrinkage"
echo "2. output.ftr = output ftr file shrinked"
echo "3. prefix     = prefix of input ftr files"
echo "4. from       = 0, 10, or 90"
echo "5. to         = 9, 19, or 99"
echo "6. sufix      = sufix"
exit 1
fi

ratio=$1; shift
output=$1; shift
prefix=$1; shift
from=$1; shift
to=$1; shift
sufix=$1; shift
range=$(../src/expand_fn.sh ${prefix} ${from} ${to} ${sufix} 0)
ds=$(../src/expand_fn.sh ${prefix} ${from} ${to} ${sufix} 1)

for f in $ds ; do
    if [ ! -e $f ]; then
    echo dataset file = $f does not exist
    exit 1
    fi
done

echo $ds

cflags="-O3 -Wall -Wno-strict-overflow"
cflags="$cflags -DCOMPILE_TIME"
cflags="$cflags -DDEEP1B"
cflags="$cflags -DFTR_ON_SECONDARY_MEMORY"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"
cflags="$cflags -DSEQUENTIAL_FILTERING"
cflags="$cflags -DPARTITION_TYPE_GHP"
cflags="$cflags -DBLOCK_SIZE=1"


gcc $cflags ../src/shrink_FTR.c ../src/ftr.c -o shrink_FTR
if [ $? == 1 ] ; then
exit 1
fi

./shrink_FTR $ratio $output $ds

