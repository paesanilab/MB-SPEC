#!/bin/sh
set -e

if [ -f Makefile ]; then
    make distclean
fi

autoreconf -fi

./configure CC=icx CXX=icpx \
    CFLAGS="-g -O2 -fopenmp -Wall -std=c99" \
    CXXFLAGS="-g -O2 -fopenmp -Wall"

make -j4
