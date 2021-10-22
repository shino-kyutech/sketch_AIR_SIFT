#!/bin/bash

s=""
for f in {00..01} ; do
s="$s base_$f.ftr"
done

echo $s

./make_answer_by_sequential_search.sh query_00.ftr 100 15 PRINTOUT $s > answer_00_01.csv


