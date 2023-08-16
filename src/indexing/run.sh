#!/bin/bash

bin/pmetindex \
    -f ../../results/PMETindex_promoters/fimo \
    -k 5 \
    -n 5000 \
    -p ../../results/PMETindex_promoters/promoter_lengths.txt \
    -o test_result