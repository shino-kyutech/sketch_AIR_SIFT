#!/bin/bash
for f in 0 1 2 3 4 ; do 
./make_sorted_ftr_on_memory.sh pivot_QBP_t1000_w28_sd00_ss10000 $(($f * 20)) $(($f * 20 + 19)) 
done
