#!/bin/bash
# without using optional Armadillo library, Optimization: DEBUG / O3-Level, no further architecture-specific optimization
#e.g. ./compile_XNDiff.sh g++ DEBUG 0 NONE NONE
#e.g. ./compile_XNDiff.sh g++ 3 0 NONE NONE
# using optional Armadillo library, Optimization: DEBUG / O3-Level, no further architecture-specific optimization
#e.g. ./compile_XNDiff.sh g++ DEBUG 1 NONE NONE
#e.g. ./compile_XNDiff.sh g++ 3 1 NONE NONE
# using optional Armadillo library, Optimization O3-Level, no further architecture-specific optimization, ignore specific gcc warnings
#e.g. ./compile_XNDiff.sh g++ 3 1 NONE "-Wno-maybe-uninitialized"
# using several -Wno does not yet work
# ./compile_XNDiff.sh g++ 3 1 NONE "-Wno-maybe-uninitialized -Wno-unused-result"
make clean
echo "make CC=$1 OLEVEL=$2 ARMADILLO=$3 OOPTIONS=$4 EXCLUDEWARN=$5"
make CC=$1 OLEVEL=$2 ARMADILLO=$3 OOPTIONS=$4 EXCLUDEWARN=$5

