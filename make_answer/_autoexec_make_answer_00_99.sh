#!/bin/bash
if [ $# -ne 1 ] ; then
echo "Usage> スクリプト.sh range"
echo "1. range    = range of dataset to be searched, 00_01, 00_09, 00_99" 
exit 1
fi

range=$1

if [ $range == 00_01 ] ; then
s=""
for f in {00..01} ; do
s="$s base_$f.ftr"
done
elif [ $range == 00_09 ] ; then
s=""
for f in {00..09} ; do
s="$s base_$f.ftr"
done
elif [ $range == 00_99 ] ; then
s=""
for f in {00..99} ; do
s="$s base_$f.ftr"
done
else
exit
fi

echo $s

exit

./make_answer_by_sequential_search.sh query_00.ftr 100 15 PRINTOUT $s > answer_${range}.csv


