#!/bin/bash

#pv=pivot_Q200_A1000000_L0_w20_s1_k20_sf0.6_100K.pivot
#pv=pivot_Q200_A100000_L0_w22_s1_k15_sf0.7_1M.pivot
pv=pivot_Q200_A100000_L0_w26_s1_k10_sf0.6_1M.pivot

for f in 0 1 2 3 4 ; do
./make_bucket.sh $pv $(($f * 20)) $(($f * 20 + 19)) 15
done

#./make_bucket.sh $pv 0 99 15

for f in 0 1 2 3 4 ; do 
./make_sorted_ftr_on_memory.sh $pv $(($f * 20)) $(($f * 20 + 19)) INTEL
done

./merge_sorted_ftr.sh $pv 0 99 20 INTEL
