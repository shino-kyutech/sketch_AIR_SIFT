#!/bin/bash

s=$RANDOM
b=800
t=8
r=1

for n in 100 200 400 800 1600 ; do

./search_by_sequential_and_random_read.sh query_00 0 99 SSD $n 100 1000 $r $b $t RANDOM $s
s=$(($s + 1))

done

exit

./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100    1  1 SEQUENTIAL 500
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100  400  1 SEQUENTIAL 501
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100  800  1 SEQUENTIAL 502
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100 1200  1 SEQUENTIAL 503

./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100    1  8 SEQUENTIAL 504
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100  400  8 SEQUENTIAL 505
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100  800  8 SEQUENTIAL 506
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100 1200  8 SEQUENTIAL 507

./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100    1 16 SEQUENTIAL 508
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100  400 16 SEQUENTIAL 509
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100  800 16 SEQUENTIAL 510
./search_by_sequential_and_random_read.sh query_00 0 99 SSD 1000000 100 100 1200 16 SEQUENTIAL 511

./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100    1  1 RANDOM     600
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100  400  1 RANDOM     601
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100  800  1 RANDOM     602
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100 1200  1 RANDOM     603

./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100    1  8 RANDOM     604
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100  400  8 RANDOM     605
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100  800  8 RANDOM     606
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100 1200  8 RANDOM     607

./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100    1 16 RANDOM     608
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100  400 16 RANDOM     609
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100  800 16 RANDOM     610
./search_by_sequential_and_random_read.sh query_00 0 99 SSD   10000 100 100 1200 16 RANDOM     611

