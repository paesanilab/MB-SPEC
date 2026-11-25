#!/bin/sh
# Build the reference commit in a dedicated checkout with its own build files.
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 HISTORICAL_CHECKOUT" >&2
    exit 2
fi
: "${MKLROOT:?initialize Intel MKL first}"

source_dir=$(cd "$(dirname "$0")" && pwd)
master=$(cd "$1" && pwd)
commit=$(git -C "$master" rev-parse HEAD)
expected=$(sed -n 's/^master commit: //p' "$source_dir/historical/references/PROVENANCE")
if [ "$commit" != "$expected" ]; then
    echo "expected historical commit $expected, found $commit" >&2
    exit 1
fi

# Clang preserves the original LAMMPS column order with the legacy C++11 flags.
cxx=${CXX:-clang++}
mkl_cflags="-I$MKLROOT/include"
mkl_libs="-L$MKLROOT/lib -Wl,-rpath,$MKLROOT/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"

echo "Building historical master ($commit) in $master"
cd "$master"
autoreconf -fi
./configure CXX="$cxx" CXXFLAGS="-g -O2 -fopenmp" \
    --with-mkl-cflags="$mkl_cflags" --with-mkl-libs="$mkl_libs"
# The old top-level Makefile expects a missing autogen.sh.
make -C src -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
