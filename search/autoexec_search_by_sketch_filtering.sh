#!/bin/bash
if [ $# -ne 16 ] ; then
echo "Usage> "
echo "1. pivot       = pivot file (without .csv)"
echo "2. range       = 00_01, 00_09, or 00_99"
echo "3. query       = 00, 01, ... , or 99"
echo "4. num_c       = number of candidates (ppm 百万分率)" 
echo "5. num_c1      = number of candidates" 
echo "6. num_c2      = number of candidates" 
echo "7. num_c3      = number of candidates" 
echo "8. num_c4      = number of candidates" 
echo "9. num_c5      = number of candidates" 
echo "10. num_k      = k (number of neighbors)"
echo "11. num_q      = number of queries (0 -> all)"
echo "12. score_p    = p of score * 10"
echo "13. #threads   = number of threads"
echo "14. mode       = 0 -> sequential filtering, 1 -> using bucket, 2 -> filtering by sketch enumeration, 3 -> filtering by sketch enumeration c2_n, 4 -> filtering by hamming"
echo "15. ftr_on     = 0 -> main memory, 1 -> secondary memory, 2 -> secondary memory (using sorted ftr), 3 -> secondary memory (using unsorted catenated ftr)"
echo "16. storage    = HDD, SSD, ADATA, INTEL"
exit 1
fi

pivot=$1; shift
range=$1; shift
query=$1; shift
nc=$1; shift
n1=$1; shift
n2=$1; shift
n3=$1; shift
n4=$1; shift
n5=$1; shift
nk=$1; shift
nq=$1; shift
sp=$1; shift
nt=$1; shift
md=$1; shift
mm=$1; shift
st=$1; shift

result="result/${pivot}_${range}.csv"

if [ $mm == 0 ] || [ $mm == 1 ]; then
	if [ $range == "00_00" ] ; then
	s="base_00.ftr"
	elif [ $range == "00_01" ] ; then
	s=""
	for f in {00..01} ; do
	s="$s base_$f.ftr"
	done
	elif [ $range == "00_09" ] ; then
	s=""
	for f in {00..09} ; do
	s="$s base_$f.ftr"
	done
	elif [ $range == "00_49" ] ; then
	s=""
	for f in {00..49} ; do
	s="$s base_$f.ftr"
	done
	elif [ $range == "00_99" ] ; then
	s=""
	for f in {00..99} ; do
	s="$s base_$f.ftr"
	done
	fi
elif [ $mm == 2 ] ; then
	s=${pivot}_${range}.sftr
else
	s=base_${range}.ftr
fi

./search_by_sketch_filtering_on_secondary_memory.sh $st $pivot $range $query $result $md $mm $nc $n1 $n2 $n3 $n4 $n5 $nk $nq $sp $nt $s
