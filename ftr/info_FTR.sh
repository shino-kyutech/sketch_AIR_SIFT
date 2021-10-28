#!/bin/bash

cflags="-O3 -Wall -Wno-strict-overflow"
cflags="$cflags -DCOMPILE_TIME"
cflags="$cflags -DDEEP1B"
cflags="$cflags -DFTR_ON_SECONDARY_MEMORY"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"
cflags="$cflags -DSEQUENTIAL_FILTERING"
cflags="$cflags -DPARTITION_TYPE_GHP"


gcc $cflags ../src/info_FTR.c ../src/ftr.c -o info_FTR
if [ $? == 1 ] ; then
exit 1
fi

./info_FTR $1

