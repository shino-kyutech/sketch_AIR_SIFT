#!/bin/bash

th=8
nq=10

#recall = 95%
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w16_sd00_ss10000_seed1 00_99 01 23320 0 0 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w18_sd00_ss10000_seed1 00_99 01  12900 0 0 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w20_sd00_ss10000_seed1 00_99 01  13500 0 0 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 01   7430 0 0 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w24_sd00_ss10000       00_99 01   7410 0 0 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w26_sd00_ss10000       00_99 01   3920 0 0 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w28_sd00_ss10000       00_99 01   3620 0 0 0 0 0 100 1000 10 $th 3 2 INTEL

#recall = 90%
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w16_sd00_ss10000_seed1 00_99 01 10700 0 23320 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w18_sd00_ss10000_seed1 00_99 01  6900 0 12900 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w20_sd00_ss10000_seed1 00_99 01  5720 0 13500 0 0 0 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 01  3610 3120 0 0 0 0 100 $nq 10 8 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 00  3610 3120 0 0 0 0 100 $nq 10 1 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 01  3610 3120 0 0 0 0 100 $nq 10 1 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 00  3610 3120 0 0 0 0 100 $nq 10 8 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w24_sd00_ss10000       00_99 01  2680 0  7410 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w26_sd00_ss10000       00_99 01  1680 0  3920 0 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w28_sd00_ss10000       00_99 01  1570 0  3620 0 0 0 100 1000 10 $th 3 2 INTEL

#recall = 85%
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w16_sd00_ss10000_seed1 00_99 01 6380 0 10700 23320 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w18_sd00_ss10000_seed1 00_99 01 3600 0  6900 12900 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w20_sd00_ss10000_seed1 00_99 01 2900 0  5720 13500 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 01 2000 0  3610  7430 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w24_sd00_ss10000       00_99 01 1540 0  2680  7410 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w26_sd00_ss10000       00_99 01 1020 0  1680  3920 0 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w28_sd00_ss10000       00_99 01  800 0  1570  3620 0 0 100 1000 10 $th 3 2 INTEL

#recall = 80%
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w16_sd00_ss10000_seed1 00_99 01 3970 0 6380 10700 23320 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w18_sd00_ss10000_seed1 00_99 01 2300 0 3600  6900 12900 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w20_sd00_ss10000_seed1 00_99 01 1710 0 2900  5720 13500 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 01 1130 0 2000  3610  7430 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w24_sd00_ss10000       00_99 01  880 0 1540  2680  7410 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w26_sd00_ss10000       00_99 01  580 0 1020  1680  3920 0 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w28_sd00_ss10000       00_99 01  490 0  800  1570  3620 0 100 1000 10 $th 3 2 INTEL

# recall = 75% (query = 01)
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w16_sd00_ss10000_seed1 00_99 01 3000 0 3970 6380 10700 23320 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w18_sd00_ss10000_seed1 00_99 01 1850 0 2300 3600  6900 12900 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w20_sd00_ss10000_seed1 00_99 01 1350 0 1710 2900  5720 13500 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 01  850 0 1130 2000  3610  7430 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w24_sd00_ss10000       00_99 01  630 0  880 1540  2680  7410 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w26_sd00_ss10000       00_99 01  430 0  580 1020  1680  3920 100 1000 10 $th 3 2 INTEL
#./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w28_sd00_ss10000       00_99 01  360 0  490  800  1570  3620 100 1000 10 $th 3 2 INTEL

exit

./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w16_sd00_ss10000_seed1 00_99 00 1610 2520 3970 6380 10700 23320 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w18_sd00_ss10000_seed1 00_99 00 1210 1630 2300 3600  6900 12900 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w20_sd00_ss10000_seed1 00_99 00  870 1250 1710 2900  5720 13500 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w22_sd00_ss10000_seed1 00_99 00  520  770 1130 2000  3610  7430 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w24_sd00_ss10000       00_99 00  400  590  880 1540  2680  7410 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w26_sd00_ss10000       00_99 00  260  370  580 1020  1680  3920 100 1000 10 $th 3 2 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_QBP_t1000_w28_sd00_ss10000       00_99 00  210  330  490  800  1570  3620 100 1000 10 $th 3 2 INTEL

./autoexec_search_by_sketch_filtering.sh pivot_PQBP_t10_w96_sd00_ss10000_np2_sr1000_seed1  00_99 38    0 0 0 0 0 100 100 10 $th 0 1 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_PQBP_t10_w192_sd00_ss10000_np8_sr1000_seed1 00_99 4.5   0 0 0 0 0 100 100 10 $th 0 1 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_PQBP_t10_w384_sd00_ss10000_np4_sr1000_seed1 00_99 0.362 0 0 0 0 0 100 100 10 $th 0 1 INTEL

./autoexec_search_by_sketch_filtering.sh pivot_PQBP_t10_w96_sd00_ss10000_np2_sr1000_seed1  00_99 4.3  7    14   38 0 0 100 100 10 $th 0 1 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_PQBP_t10_w192_sd00_ss10000_np8_sr1000_seed1 00_99 0.35 0.65 1.2 4.5 0 0 100 100 10 $th 0 1 INTEL
./autoexec_search_by_sketch_filtering.sh pivot_PQBP_t10_w384_sd00_ss10000_np4_sr1000_seed1 00_99 0.043 0.073 0.134 0.362 0 0 100 100 10 $th 0 1 INTEL
