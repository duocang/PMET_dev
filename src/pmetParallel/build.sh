#!/bin/bash

# 创建目标目录
mkdir -p bin

cd bin
# 执行CMake命令生成Makefile
cmake -DCMAKE_BUILD_TYPE=Debug ..

make

sleep 1
rm Makefile
rm -rf CMake*
rm -rf cmake_install.cmake

cd ..

# remove genes not present in pre-computed pmet index (universe.txt)
grep -Ff test_data/universe.txt test_data/gene.txt > test_data/gene.txttemp

# Run PMET
bin/pmetParallel \
    -d . \
    -g test_data/gene.txttemp \
    -i 4 \
    -p test_data/promoter_lengths.txt \
    -b test_data/binomial_thresholds.txt \
    -c test_data/IC.txt \
    -f test_data/fimohits \
    -o test_result \
    -t 2

rm test_data/gene.txttemp

cat test_result/*.txt > test_result/motif_output.txt
rm test_result/temp*.txt
