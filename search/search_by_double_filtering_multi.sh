#!/bin/bash
if [ $# -ne 19 ] ; then
echo "Usage>"
echo "1.  storage    = HDD, SSD, ADATA, INTEL"
echo "2.  1st pivot  = pivot file for narrow sketch (for 1st filtering)"
echo "3.  2nd pivot  = pivot file for expanded sketch (for 2nd filtering)"
echo "4.  range      = range of ftr files" 
echo "5.  1st query  = file number of 1st query file" 
echo "6.  2nd query  = file number of 2nd query file (-1 for skip)" 
echo "7.  3rd query  = file number of 3rd query file (-1 for skip)" 
echo "8.  result.csv = result file" 
echo "9.  nc1_start  = start: number of candidates of 1st filtering (ppm)"
echo "10. nc1_end    = end "
echo "11. nc1_step   = step"
echo "12. nc2_start  = start: number of candidates of 2nd filtering (ppm)"
echo "13. nc2_end    = end"
echo "14. nc2_step   = step"
echo "15. num_k      = k (number of neighbors)"
echo "16. num_q      = number of queries (0 -> all)"
echo "17. p_1st      = p of score * 10 for 1st filtering"
echo "18. p_2nd      = p of score * 10 for 2nd filtering"
echo "19. #threads   = number of threads"

#echo "6.  mode       = 0 -> sequential filtering, 1 -> using bucket, 2 -> filtering by sketch enumeration, 3 -> filtering by sketch enumeration c2_n, 4 -> filtering by hamming"
#echo "7.  ftr_on     = 0 -> main memory, 1 -> secondary memory, 2 -> secondary memory (using sorted ftr), 3 -> secondary memory (using unsorted catenated ftr)"
#echo "18 - dataset_1.ftr, ... , dataset_n file = dataset files" 

exit 1
fi

#set -x

ulimit -Sn 10000
#ulimit -a

st=$1; shift
if [ $st == HDD ] || [ $st == SSD ] || [ $st == ADATA ] || [ $st == INTEL ] ; then
source ../src/set_directory.sh $st
else
echo "invalid storage $st: storage = HDD, SSD, ADATA, or INTEL"
exit
fi

pivot1=$1; shift

p1="$pv_dir/${pivot1}.csv"
if [ ! -e $p1 ]; then
  echo pivot file for 1st filtering = $p1 does not exist.
  exit
fi

pivot2=$1; shift

p2="$pv_dir/${pivot2}.csv"
if [ ! -e $p2 ]; then
  echo pivot file for 2nd filtering = $p2 does not exist.
  exit
fi

range=$1; shift
b1="$bk_dir/${pivot1}_${range}.bkt"
if [ ! -e $b1 ]; then
  echo bucket file for 1st filtering = $b1 does not exist.
  exit
fi

b2="$bk_dir/${pivot2}_${range}.bkt"
if [ ! -e $b2 ]; then
  echo bucket file for 2nd filtering = $b2 does not exist.
  exit
fi

w1=$(../src/pivot_property.sh -w $p1)
echo NARROW_PJT_DIM = $w1
dm=$(../src/pivot_property.sh -d $p1)
echo FTR_DIM = $dm
pt1=$(../src/pivot_property.sh -p $p1)
echo partition type of 1st pivot = $pt1

w2=$(../src/pivot_property.sh -w $p2)
echo EXPANDED_PJT_DIM = $w2
d2=$(../src/pivot_property.sh -d $p2)
pt=$(../src/pivot_property.sh -p $p2)
np=$(../src/pivot_property.sh -n $p2)
echo partition type of 2nd pivot = $pt, number of partitioned spaces = $np

if [ $dm -ne $d2 ] ; then
  echo ftr dimensions of 1st and 2nd pivots NOT EQUAL.
  exit 
fi

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

q2nd=$1; shift
if [ $q2nd -eq -1 ] ; then
  q2=NONE
  a2=NONE
else
  q2="$qr_dir/query_${q2nd}.ftr"
  if [ ! -e $q2 ]; then
    echo query file = $q2 does not exist.
    exit
  fi
  a2="$qr_dir/answer_q_${q2nd}_r_${range}.csv"
  if [ ! -e $a2 ]; then
    echo answer file = $a2 does not exist.
    exit
  fi
fi

q3rd=$1; shift
if [ $q3rd -eq -1 ] ; then
  q3=NONE
  a3=NONE
else
  q3="$qr_dir/query_${q3rd}.ftr"
  if [ ! -e $q3 ]; then
    echo query file = $q3 does not exist.
    exit
  fi
  a3="$qr_dir/answer_q_${q3rd}_r_${range}.csv"
  if [ ! -e $a3 ]; then
    echo answer file = $a3 does not exist.
    exit
  fi
fi

rs="$1"; shift

ft="DOUBLE_FILTERING"
fo="SECONDARY_MEMORY"
fr="ARRANGEMENT_ASIS"

nc1_start="$1"; shift
nc1_end="$1"; shift
nc1_step="$1"; shift
nc2_start="$1"; shift
nc2_end="$1"; shift
nc2_step="$1"; shift
nk="$1"; shift
nq="$1"; shift
sp1="$1"; shift
sp2="$1"; shift
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
cflags0="$cflags0 -DNARROW_PJT_DIM=$w1"
cflags0="$cflags0 -DEXPANDED_PJT_DIM=$w2"
#cflags0="$cflags0 -DEXPANDED_TABLE_WIDE"

cflags0="$cflags0 -DPRIORITY=1"

cflags0="$cflags0 -DSCORE_P_1ST=${sp1}"
cflags0="$cflags0 -DSCORE_P_2ND=${sp2}"

#cflags0="$cflags0 -DPARTITION_TYPE_QBP"

#if [ $pt == 3 ] ; then
cflags0="$cflags0 -DPARTITION_TYPE_QBP_AND_PQBP"
cflags0="$cflags0 -DNUM_PART=$np"
#cflags0="$cflags0 -DNUM_PART=1"
#else
#fi

cflags0="$cflags0 -DFTR_ON_${fo}"
cflags0="$cflags0 -DFTR_${fr}"

cflags="$cflags -DBLOCK_SIZE=800"
cflags="$cflags -D${ft}"

cflags="$cflags -DNARROW_PIVOT_FILE=\"$p1\""
cflags="$cflags -DNARROW_BUCKET_FILE=\"$b1\""
cflags="$cflags -DEXPANDED_PIVOT_FILE=\"$p2\""
cflags="$cflags -DEXPANDED_BUCKET_FILE=\"$b2\""
cflags="$cflags -DQUERY_FILE=\"$qr\""
cflags="$cflags -DANSWER_FILE=\"$an\""
cflags="$cflags -DRESULT_FILE=\"result/$rs\""
cflags="$cflags -DQUERY_2ND_FILE=\"$q2\""
cflags="$cflags -DANSWER_2ND_FILE=\"$a2\""
cflags="$cflags -DQUERY_3RD_FILE=\"$q3\""
cflags="$cflags -DANSWER_3RD_FILE=\"$a3\""

cflags="$cflags -DNUM_CANDIDATES_1ST=$nc1_start"
cflags="$cflags -DNUM_CANDIDATES_1ST_END=$nc1_end"
cflags="$cflags -DNUM_CANDIDATES_1ST_STEP=$nc1_step"

cflags="$cflags -DNUM_CANDIDATES_2ND=$nc2_start"
cflags="$cflags -DNUM_CANDIDATES_2ND_END=$nc2_end"
cflags="$cflags -DNUM_CANDIDATES_2ND_STEP=$nc2_step"

cflags="$cflags -DNUM_K=$nk"

#cflags="$cflags -DSELECT_K_BY_QUICK_SELECT_K"
cflags="$cflags -DSELECT_K_BY_QUICK_SELECT_K_PARA_WORK"
#cflags="$cflags -DSELECT_K_BY_QUICK_SORT"


if [ $nq > 0 ] ; then
cflags="$cflags -DNUM_Q=$nq"
fi

echo $cflags0
echo $cflags

ds=$(../src/expand_filenames.sh base_ $range .ftr)

files=""
for f in $ds ; do
	if [ ! -e $ds_dir/$f ]; then
	    echo dataset file = $ds_dir/$f does not exist
	    exit 1
	fi
	files="$files $ds_dir/$f"
done

#echo $files

l1=bit_op
l2=ftr
l3=kNN_search
l4=double_sketch
l5=quick
l6=e_time

lp="$pr_dir/$l1.c $pr_dir/$l2.c $pr_dir/$l3.c $pr_dir/$l4.c $pr_dir/$l5.c $pr_dir/$l6.c"

echo gcc $cflags0 $cflags -c $lp
gcc $cflags0 $cflags -c $lp

if [ $? == 1 ] ; then
echo compile error for libraly
exit 1
fi

lb="$l1.o $l2.o $l3.o $l4.o $l5.o $l6.o -lm"

pr=search_by_double_filtering_multi

gcc $cflags0 $cflags $pr_dir/$pr.c $lb -o $pr

if [ $? == 1 ] ; then
exit 1
fi

time ./$pr $files

