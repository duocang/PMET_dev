#!/bin/bash

mkdir -p test_result

# remove genes not present in pre-computed pmet index (universe.txt)
grep -Ff ../../results/01_homotypic_promoters/universe.txt ../../data/gene.txt > test_result/gene.txttemp

# Run PMET
bin/pmetParallel \
    -d . \
    -g test_result/gene.txttemp \
    -i 4 \
    -p ../../results/01_homotypic_promoters/promoter_lengths.txt  \
    -b ../../results/01_homotypic_promoters/binomial_thresholds.txt  \
    -c ../../results/01_homotypic_promoters/IC.txt  \
    -f ../../results/01_homotypic_promoters/fimohits  \
    -t 8 \
    -o test_result

rm test_result/gene.txttemp

cat test_result/*.txt > test_result/motif_output.txt
rm test_result/temp*.txt
