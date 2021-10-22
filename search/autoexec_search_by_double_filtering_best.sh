#!/bin/bash

# RAM = 20GB (process = 15GB)

p2=pivot_PQBP_t10_w192_sd00_ss10000_np8_sr1000_seed1
range=00_99
query=02
nk=100
nq=1000
s1=10
s2=10
nt=8

p1=pivot_QBP_t1000_w22_sd00_ss10000_seed1
result=double_filtering_best_w22_q00.csv

./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 2500 2500 1000 0.6 0.6 0.2 $nk $nq $s1 $s2 $nt
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 8000 8000 1000 1.6 1.6 0.1 $nk $nq $s1 $s2 $nt

p1=pivot_QBP_t1000_w24_sd00_ss10000
result=double_filtering_best_w24_q00.csv

./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 2000 2000 1000 0.6 0.6 0.2 $nk $nq $s1 $s2 $nt
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 8000 8000 1000 1.6 1.6 0.1 $nk $nq $s1 $s2 $nt
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 7500 7500 1000 1.6 1.6 0.1 $nk $nq $s1 $s2 $nt

exit

p1=pivot_QBP_t1000_w26_sd00_ss10000
result=double_filtering_best_w26_q00.csv

./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 2000 2000 1000 0.4 0.4 0.2 $nk $nq $s1 $s2 $nt
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 4000 4000 1000 1.6 1.6 0.1 $nk $nq $s1 $s2 $nt

p1=pivot_QBP_t1000_w28_sd00_ss10000
result=double_filtering_best_w28_q00.csv

./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 2000 2000 1000 0.4 0.4 0.2 $nk $nq $s1 $s2 $nt
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 3000 3000 1000 1.6 1.6 0.2 $nk $nq $s1 $s2 $nt

