#!/bin/bash
#without using optional Armadillo library, Optimization: DEBUG / O3-Level, no further architecture-specific optimization
#e.g. ./compile_XNDiff.sh g++ DEBUG 0 NONE
#e.g. ./compile_XNDiff.sh g++ 3 0 NONE
#using optional Armadillo library, Optimization: DEBUG / O3-Level, no further architecture-specific optimization
#e.g. ./compile_XNDiff.sh g++ DEBUG 1 NONE
#e.g. ./compile_XNDiff.sh g++ 3 1 NONE
make clean
echo "make CC=$1 OLEVEL=$2 ARMADILLO=$3 OOPTIONS=$4"
make CC=$1 OLEVEL=$2 ARMADILLO=$3 OOPTIONS=$4

