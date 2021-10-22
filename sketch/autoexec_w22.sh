#!/bin/bash

w=18
tr=1000
sd=00
ss=10000
seed=1
./make_sketch_by_QBP.sh $w $tr $sd $ss 15 $seed

pv="pivot_QBP_t${tr}_w${w}_sd${sd}_ss${ss}_seed${seed}"

for f in 0 1 2 3 4 ; do
./make_bucket.sh ${pv} $(($f * 20)) $(($f * 20 + 19)) 15
done

./make_bucket.sh ${pv} 0 99 15

./num_k_for_recall.sh ${pv} 0.4 1.5 0.1 00_99 00 15 > score_p_w${w}_00_99.csv

for f in 0 1 2 3 4 ; do 
./make_sorted_ftr_on_memory.sh ${pv} $(($f * 20)) $(($f * 20 + 19)) 
done

./merge_sorted_ftr.sh ${pv} 0 99 20
