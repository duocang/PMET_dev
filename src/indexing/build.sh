#!/bin/bash

mkdir -p bin

cd bin

cmake -DCMAKE_BUILD_TYPE=Debug ..

make

# sleep 1
# rm Makefile
# rm -rf CMake*
# rm -rf cmake_install.cmake

cd ..
mkdir -p test_result/fimohits
