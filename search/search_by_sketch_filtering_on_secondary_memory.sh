#!/bin/bash
if [ $# -lt 17 ] ; then
echo "Usage>"
echo "1.  storage    = HDD, SSD, ADATA, INTEL"
echo "2.  pivot      = pivot file"
echo "3.  range      = range of ftr files" 
echo "4.  query      = query file" 
echo "5.  result.csv = result file" 
echo "6.  mode       = 0 -> sequential filtering, 1 -> using bucket, 2 -> filtering by sketch enumeration, 3 -> filtering by sketch enumeration c2_n, 4 -> filtering by hamming"
echo "7.  ftr_on     = 0 -> main memory, 1 -> secondary memory, 2 -> secondary memory (using sorted ftr), 3 -> secondary memory (using unsorted catenated ftr)"
echo "8.  num_can    = number of candidates (ppm)"
echo "9.  num_can1   = number of candidates (ppm), 0 for end"
echo "10. num_can2   = number of candidates (ppm)"
echo "11. num_can3   = number of candidates (ppm)"
echo "12. num_can4   = number of candidates (ppm)"
echo "13. num_can5   = number of candidates (ppm)"
echo "14. num_k      = k (number of neighbors)"
echo "15. num_q      = number of queries (0 -> all)"
echo "16. score_p    = p of score * 10"
echo "17. #threads   = number of threads"
echo "18 - dataset_1.ftr, ... , dataset_n file = dataset files" 
exit 1
fi

#set -x

ulimit -Sn 4000

st=$1; shift
if [ $st == HDD ] || [ $st == SSD ] || [ $st == ADATA ] || [ $st == INTEL ] ; then
source ../src/set_directory.sh $st
else
echo "invalid storage $st: storage = HDD, SSD, ADATA, or INTEL"
exit
fi

pivot=$1; shift

pv="$pv_dir/${pivot}.csv"
if [ ! -e $pv ]; then
  echo pivot file = $pv does not exist.
  exit
fi

range=$1; shift
bk="$bk_dir/${pivot}_${range}.bkt"
if [ ! -e $bk ]; then
  echo bucket file = $bk does not exist.
  exit
fi
wd=$(${pr_dir}/pivot_property.sh -w $pv)
dm=$(${pr_dir}/pivot_property.sh -d $pv)
pt=$(${pr_dir}/pivot_property.sh -p $pv)
np=$(${pr_dir}/pivot_property.sh -n $pv)

query=$1; shift
qr="$qr_dir/query_${query}.ftr"
if [ ! -e $qr ]; then
  echo query file = $qr does not exist.
  exit
fi

an="$qr_dir/answer_q_${query}_r_${range}.csv"
if [ ! -e $an ]; then
  echo answer file = $an does not exist.
  exit
fi

rs="$1"; shift
md="$1"; shift

if [ $md == 0 ] ; then
ft="SEQUENTIAL_FILTERING"
elif [ $md == 1 ] ; then
ft="SEQUENTIAL_FILTERING_USING_BUCKET"
elif [ $md == 2 ] ; then
ft="FILTERING_BY_SKETCH_ENUMERATION"
elif [ $md == 3 ] ; then
ft="FILTERING_BY_SKETCH_ENUMERATION_C2N"
elif [ $md == 4 ] ; then
ft="SEQUENTIAL_FILTERING_USING_HAMMING"
fi

mm="$1"; shift

if [ $mm == 0 ] ; then
	fo="MAIN_MEMORY"
	fr="ARRANGEMENT_ASIS"
elif [ $mm == 1 ] ; then
	fo="SECONDARY_MEMORY"
	fr="ARRANGEMENT_ASIS"
elif [ $mm == 2 ] ; then
	fo="SECONDARY_MEMORY"
	fr="SORT_BY_SKETCH_IN_ADVANCE"
elif [ $mm == 3 ] ; then
	fo="SECONDARY_MEMORY"
	fr="ARRANGEMENT_ASIS"
fi

nc="$1"; shift
n1="$1"; shift
n2="$1"; shift
n3="$1"; shift
n4="$1"; shift
n5="$1"; shift
nk="$1"; shift
nq="$1"; shift
sp="$1"; shift
nt="$1"; shift

#if [ $nt == 1 ] ; then
#cflags0="-O3 -Wall -Wno-strict-overflow"
#else
cflags0="-O3 -fopenmp -Wall -Wno-strict-overflow"
cflags0="$cflags0 -DNUM_THREADS=$nt"
#fi

cflags0="$cflags0 -DCOMPILE_TIME"
cflags0="$cflags0 -DDEEP1B"
cflags0="$cflags0 -DFTR_DIM=$dm"
cflags0="$cflags0 -DPJT_DIM=$wd"

if [ $wd -gt 64 ] ; then
cflags0="$cflags0 -DEXPANDED_SKETCH"
elif [ $wd -gt 32 ] ; then
cflags0="$cflags0 -DWIDE_SKETCH"
else
cflags0="$cflags0 -DNARROW_SKETCH"
fi

cflags0="$cflags0 -DPRIORITY=1"

if [ $sp == 10 ] ; then
cflags0="$cflags0 -DSCORE_1"
elif [ $sp == 20 ] ; then
cflags0="$cflags0 -DSCORE_2"
else
cflags0="$cflags0 -DSCORE_P=${sp}"
fi

if [ $pt == 3 ] ; then
cflags0="$cflags0 -DPARTITION_TYPE_PQBP"
cflags0="$cflags0 -DNUM_PART=$np"
else
cflags0="$cflags0 -DPARTITION_TYPE_QBP"
fi

cflags0="$cflags0 -DFTR_ON_${fo}"
cflags0="$cflags0 -DFTR_${fr}"

cflags="$cflags -DBLOCK_SIZE=800"
#cflags="$cflags -DFTR_SORT_BY_SKETCH"
#cflags="$cflags -DFTR_SORT_BY_SKETCH_IN_ADVANCE"
#cflags="$cflags -DSEQUENTIAL_FILTERING"
#cflags="$cflags -DSEQUENTIAL_FILTERING_USING_BUCKET"
cflags="$cflags -D${ft}"
#cflags="$cflags -DFILTERING_BY_SKETCH_ENUMERATION"

cflags="$cflags -DPIVOT_FILE=\"$pv\""
cflags="$cflags -DBUCKET_FILE=\"$bk\""
cflags="$cflags -DQUERY_FILE=\"$qr\""
cflags="$cflags -DANSWER_FILE=\"$an\""
cflags="$cflags -DRESULT_FILE=\"$rs\""
cflags="$cflags -DNUM_CANDIDATES=$nc"
cflags="$cflags -DNUM_CANDIDATES1=$n1"
cflags="$cflags -DNUM_CANDIDATES2=$n2"
cflags="$cflags -DNUM_CANDIDATES3=$n3"
cflags="$cflags -DNUM_CANDIDATES4=$n4"
cflags="$cflags -DNUM_CANDIDATES5=$n5"
cflags="$cflags -DNUM_K=$nk"

if [ $nq > 0 ] ; then
cflags="$cflags -DNUM_Q=$nq"
fi

echo $cflags0
echo $cflags

if [ $mm == 0 ] || [ $mm == 1 ] || [ $mm == 3 ] ; then
	ds=" "
	count=1
	while [ "$#" -ge "1" ]; do
	  if [ ! -e $ds_dir/$1 ]; then
	    echo dataset file = $ds_dir/$1 does not exist
	    exit 1
	  fi
	  ds="$ds $ds_dir/$1 "
	  shift
	  let count=$count+1
	done
else
	ds=" "
	count=1
	while [ "$#" -ge "1" ]; do
		if [ $count -gt 1 ] ; then
			echo "too many dataset files: only one sorted dataset file should be used"
			exit 1
		fi
		if [ ! -e $sf_dir/$1 ]; then
			echo "sorted dataset file = $sf_dir/$1 does not exist"
			exit 1
		fi
		ds="$ds $sf_dir/$1 "
		shift
		let count=$count+1
	done
fi

echo $ds

l1=bit_op
l2=ftr
l3=kNN_search
l4=sketch
l5=quick
l6=e_time

lp="$pr_dir/$l1.c $pr_dir/$l2.c $pr_dir/$l3.c $pr_dir/$l4.c $pr_dir/$l5.c $pr_dir/$l6.c"

gcc $cflags0 $cflags -c $lp

if [ $? == 1 ] ; then
echo compile error for libraly
exit 1
fi

lb="$l1.o $l2.o $l3.o $l4.o $l5.o $l6.o -lm"

pr=search_by_sketch_filtering_on_secondary_memory

gcc $cflags0 $cflags $pr_dir/$pr.c $lb -o $pr
#gcc $cflags0 $cflags $pr_dir/$pr.c $lp -o $pr

if [ $? == 1 ] ; then
exit 1
fi

time ./$pr $ds

