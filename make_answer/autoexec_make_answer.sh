#!/bin/bash
if [ $# -ne 2 ] ; then
echo "Usage> "
echo "1. query    = query file number, 00, 01, ... 09" 
echo "2. range    = range of dataset to be searched, 00_01, 00_09, 00_99" 
exit 1
fi

query=$1; shift
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
elif [ $range == 00_49 ] ; then
s=""
for f in {00..49} ; do
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

./make_answer_by_sequential_search.sh query_${query}.ftr 100 15 PRINTOUT $s > answer_q_${query}_r_${range}.csv


