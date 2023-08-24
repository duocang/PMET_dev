#!/bin/bash

mkdir -p test_result

# remove genes not present in pre-computed pmet index (universe.txt)
grep -Ff ../../data/test_for_pmet/universe.txt ../../data/test_for_pmet/gene.txt > test_result/gene.txttemp

# Run PMET
bin/pmetParallel \
    -d . \
    -g test_result/gene.txttemp \
    -i 4 \
    -p ../../data/test_for_pmet/promoter_lengths.txt  \
    -b ../../data/test_for_pmet/binomial_thresholds.txt  \
    -c ../../data/test_for_pmet/IC.txt  \
    -f ../../data/test_for_pmet/fimohits  \
    -t 2 \
    -o test_result

rm test_result/gene.txttemp

cat test_result/*.txt > test_result/motif_output.txt
rm test_result/temp*.txt
