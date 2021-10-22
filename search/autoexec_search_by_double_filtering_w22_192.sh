#!/bin/bash

# RAM = 20GB (process = 15GB)

p1=pivot_QBP_t1000_w22_sd00_ss10000_seed1
#p1=pivot_QBP_t1000_w28_sd00_ss10000
p2=pivot_PQBP_t10_w192_sd00_ss10000_np8_sr1000_seed1
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
#for q in 00 50 90 ; do
#query=query_$q
./search_by_double_filtering_multi.sh INTEL $p1 $p2 $range $query $result 3000 18000 3000 2 10 2 $nk $nq $s1 $s2 $nt
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
