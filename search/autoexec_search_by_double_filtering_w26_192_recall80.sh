#!/bin/bash

# RAM = 20GB (process = 15GB)

p1=pivot_QBP_t1000_w26_sd00_ss10000
#p1=pivot_QBP_t1000_w22_sd00_ss10000_seed1
#p1=pivot_QBP_t1000_w28_sd00_ss10000
p2=pivot_PQBP_t10_w192_sd00_ss10000_np8_sr1000_seed1
range=00_99
query=00
result=double_filtering_w26_w192_q00.csv
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
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 1000 1000 200 0.8 0.8 0.05 $nk $nq $s1 $s2 $nt
#./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 1500 1500 200 0.5 0.5 0.05 $nk $nq $s1 $s2 $nt
#./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 2000 2000 250 0.4 0.4 0.05 $nk $nq $s1 $s2 $nt
#./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 5000 5000 250 0.3 0.3 0.05 $nk $nq $s1 $s2 $nt
#./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range 00 01 02 $result 3500 5000 500 0.15 0.25 0.05 $nk $nq $s1 $s2 $nt
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
