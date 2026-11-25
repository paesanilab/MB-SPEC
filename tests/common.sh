# Shared helpers for the workflow tests. Sourced by the *.test scripts.
#
# Environment:
#   srcdir                 directory containing this file (set by Automake)
#   MBSPEC_BUILDDIR        top build directory holding the built programs
#                          (default: $top_builddir, set by tests/Makefile.am)
#   MBSPEC_PROGRAM_NAMES   "current" (mbspec-<workflow>-<stage>, default) or
#                          "legacy" (<stage>, as built by master)
#   MBSPEC_GENERATE_DIR    if set, save outputs there instead of comparing
#                          them (used by generate-reference.sh)
#   MBSPEC_KEEP_WORK       if set, keep the temporary work directory
set -eu

: "${srcdir:=$(dirname "$0")}"
srcdir=$(cd "$srcdir" && pwd)
MBSPEC_BUILDDIR=$(cd "${MBSPEC_BUILDDIR:-${top_builddir:-$srcdir/..}}" && pwd)
MBSPEC_PROGRAM_NAMES=${MBSPEC_PROGRAM_NAMES:-current}
export OMP_NUM_THREADS=1

data=$srcdir/data
workflow=$(basename "$0" .test)
reference=$srcdir/reference/$workflow
work=$(mktemp -d "${TMPDIR:-/tmp}/mbspec-$workflow.XXXXXX")
failures=0

cleanup() {
    if [ -n "${MBSPEC_KEEP_WORK:-}" ]; then
        echo "work directory kept: $work"
    else
        rm -rf "$work"
    fi
}
trap cleanup EXIT
cd "$work"

# Path of a built program: prog <workflow> <stage>
prog() {
    case $1 in
        ir) dir=infrared ;;
        *) dir=$1 ;;
    esac
    if [ "$MBSPEC_PROGRAM_NAMES" = legacy ]; then
        path=$MBSPEC_BUILDDIR/src/$dir/$2
    else
        path=$MBSPEC_BUILDDIR/src/$dir/mbspec-$1-$2
    fi
    if [ ! -x "$path" ]; then
        echo "missing program: $path" >&2
        exit 99
    fi
    printf '%s\n' "$path"
}

# Run a stage, writing stdout to OUTPUT and stderr to OUTPUT.err:
#   run OUTPUT program args...
run() {
    output=$1
    shift
    echo "+ $* > $output"
    if ! "$@" > "$output" 2> "$output.err"; then
        echo "FAILED: $* (stderr follows)" >&2
        cat "$output.err" >&2
        exit 1
    fi
}

# Check or save one output file, with optional relative and absolute tolerances.
check_file() {
    name=$1
    shift
    if [ -n "${MBSPEC_GENERATE_DIR:-}" ]; then
        mkdir -p "$MBSPEC_GENERATE_DIR/$workflow"
        cp "$name" "$MBSPEC_GENERATE_DIR/$workflow/$name"
    elif ! sh "$srcdir/compare.sh" "$reference/$name" "$name" "$@"; then
        failures=$((failures + 1))
    fi
}

check() {
    for file in "$@"; do
        check_file "$file"
    done
}

finish() {
    if [ "$failures" -ne 0 ]; then
        echo "$workflow: $failures output files differ from the master reference" >&2
        exit 1
    fi
    exit 0
}
