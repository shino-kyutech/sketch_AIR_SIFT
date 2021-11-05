#!/bin/bash

# RAM = 20GB (process = 15GB)

p1=pivot_Q1000_A100000_L0_w22_s1_k20_1M.pivot
#p1=pivot_QBP_t1000_w22_sd00_ss10000_seed1
#p1=pivot_QBP_t1000_w28_sd00_ss10000
p2=pivot_FS_Q10_A10000_L0_w192_n1_s1_k20_100K.pivot
#p2=pivot_PQBP_t10_w192_sd00_ss10000_np8_sr1000_seed1
range=00_99
query=00
result=double_filtering_w22_w192_q00.csv
k1=10000
#k2=100
nk=100
nq=1000
s1=10
s2=10
nt=8

#for nt in 8 16 ; do
#for q in 01 02 ; do
#query=query_$q
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range $query 01 02 $result 2500 3000 100 0.7 1.2 0.1 $nk $nq $s1 $s2 $nt
#./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range $query 01 02 $result 2100 2400 200 0.7 1.0 0.1 $nk $nq $s1 $s2 $nt
#done
#done

exit

for nt in 8 16 ; do
for q in 00 50 90 ; do
query=fc6_q_$q
./search_by_double_filtering_multi.sh INTEL $p1 $p2 00_96 $query $result 25000 25000 1000 3 3 1 $nk $nq $s1 $s2 $nt
done
done

for nt in 8 16 ; do
for q in 00 50 90 ; do
query=fc6_q_$q
./search_by_double_filtering_multi.sh INTEL $p1 $p2 00_96 $query $result 50000 50000 1000 6 6 1 $nk $nq $s1 $s2 $nt
done
done
