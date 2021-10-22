#!/bin/bash


range=00_09
for w in 384 768 1536 ; do
for n in 1 2 4 6 8 12 16 ; do
./make_sketch_by_PQBP.sh ${w} $n 10 00 10000 1000 15 1
pv=pivot_PQBP_t10_w${w}_sd00_ss10000_np${n}_sr1000_seed1
./make_bucket.sh ${pv} ${range} 15
./num_k_for_recall_by_SF.sh ${pv} 0.8 2.0 0.2 00_09 00 15 > score_p_w${w}_np${n}_${range}.csv 
done
done
