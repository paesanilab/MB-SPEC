#!/bin/sh
# Build an unmodified master checkout and record its workflow-test outputs as
# the reference data in tests/reference/.
#
# Usage: tests/generate-reference.sh MASTER_SOURCE OUTPUT_DIR [MASTER_COMMIT]
#
#   MASTER_SOURCE  an unmodified checkout or export of the master branch; it is
#                  built in place
#   OUTPUT_DIR     where to write <workflow>/<file>.dat and PROVENANCE
#   MASTER_COMMIT  commit id to record; defaults to `git rev-parse HEAD` in
#                  MASTER_SOURCE
#
# Master supports only Intel MKL. Initialize MKL first (e.g. source
# /opt/intel/oneapi/setvars.sh) so MKLROOT is set.
#
# Master's build files are used as they are, with two invocation details:
#   * its root Makefile.am lists a nonexistent autogen.sh, so only src/ is
#     built (`make -C src`);
#   * it never adds OpenMP flags itself (its compile.sh passed -fopenmp in
#     CXXFLAGS), so -fopenmp is passed in CXXFLAGS here too.
# No master source file is changed.
#
# The compiler defaults to clang++. Master builds src/tools with -std=c++11,
# where read-lammpstrj.cpp's `columns[column] = columns.size();` is
# unsequenced: GCC numbers LAMMPS atom columns from 1, which reads the wrong
# columns (or rejects headers whose last column is required). Clang, and
# master's documented compiler icpx (LLVM-based), number them from 0 as
# intended. C++17, used by the modern build, makes the order well defined.
set -eu

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "usage: $0 MASTER_SOURCE OUTPUT_DIR [MASTER_COMMIT]" >&2
    exit 2
fi
if [ -z "${MKLROOT:-}" ]; then
    echo "$0: MKLROOT is not set; initialize Intel MKL first" >&2
    exit 1
fi

tests=$(cd "$(dirname "$0")" && pwd)
master=$(cd "$1" && pwd)
mkdir -p "$2"
out=$(cd "$2" && pwd)
commit=${3:-$(git -C "$master" rev-parse HEAD)}
cxx=${CXX:-clang++}
cxxflags="-g -O2 -fopenmp"
mkl_cflags="-I$MKLROOT/include"
mkl_libs="-L$MKLROOT/lib -Wl,-rpath,$MKLROOT/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"

echo "Building master ($commit) in $master"
(
    cd "$master"
    autoreconf -fi
    ./configure CXX="$cxx" CXXFLAGS="$cxxflags" \
        --with-mkl-cflags="$mkl_cflags" --with-mkl-libs="$mkl_libs"
    make -C src -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
)

rm -rf "${out:?}"/ir "$out"/raman "$out"/sfg "$out"/power
for workflow in ir raman sfg power; do
    echo "Running $workflow workflow"
    srcdir=$tests MBSPEC_BUILDDIR=$master MBSPEC_PROGRAM_NAMES=legacy \
        MBSPEC_GENERATE_DIR=$out sh "$tests/$workflow.test"
done

mkl_version=$(sed -n 's/^#define INTEL_MKL_VERSION *//p' "$MKLROOT/include/mkl_version.h" 2>/dev/null || true)
cat > "$out/PROVENANCE" <<EOF
Reference outputs generated from unmodified master by tests/generate-reference.sh.
master commit: $commit
compiler: $("$cxx" --version | head -n 1)
CXXFLAGS: $cxxflags
MKL: ${mkl_version:-unknown} ($MKLROOT)
MKL libraries: $mkl_libs
OMP_NUM_THREADS: 1
EOF
echo "Reference outputs written to $out"
