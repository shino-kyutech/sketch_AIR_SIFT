#!/bin/bash

for f in 0 1 2 3 4 ; do
./make_bucket.sh $1 $(($f * 20)) $(($f * 20 + 19)) 15
done

#./make_bucket.sh $1 0 99 15

for f in 0 1 2 3 4 ; do 
./make_sorted_ftr_on_memory.sh $1 $(($f * 20)) $(($f * 20 + 19)) 
done

./merge_sorted_ftr.sh $1 0 99 20
