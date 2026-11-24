#!/bin/sh
set -eu

cd "$(dirname "$0")"

# Rebuild when compilers or flags change.
if [ -f Makefile ]; then
    make distclean
fi

./autogen.sh
./configure CXX="${CXX:-icpx}" CXXFLAGS="${CXXFLAGS--O2}" "$@"
make -j"${JOBS:-4}"
