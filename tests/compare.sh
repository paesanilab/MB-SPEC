#!/bin/sh
# Compare a program output file against a reference file.
#
# Usage: compare.sh REFERENCE OUTPUT [RTOL] [ATOL]
#
# Both files must have the same number of lines and the same number of
# whitespace-separated fields per line. Non-numeric fields must match exactly.
# Numeric fields a (reference) and b (output) match when
#     |a - b| <= RTOL * max(|a|, |b|) + ATOL * scale(column) + ulp
# where ulp is one unit in the last printed digit of the coarser of the two
# tokens (programs print 6-10 significant digits, so a last-bit difference
# can flip the final printed digit), and scale(column) is the largest
# magnitude in that column of the
# reference (or of the whole reference file, for an all-zero column), so
# near-zero values are judged relative to the column's size.
# Defaults: RTOL = 1e-8, ATOL = 1e-10. The reported maximum deviation is
# |a - b| / scale(column), the largest difference relative to its column.
set -eu

if [ $# -lt 2 ] || [ $# -gt 4 ]; then
    echo "usage: $0 REFERENCE OUTPUT [RTOL] [ATOL]" >&2
    exit 2
fi

reference=$1
output=$2
rtol=${3:-1e-8}
atol=${4:-1e-10}

for file in "$reference" "$output"; do
    if [ ! -f "$file" ]; then
        echo "compare.sh: missing file: $file" >&2
        exit 1
    fi
done
if [ ! -s "$reference" ]; then
    echo "compare.sh: empty reference file: $reference" >&2
    exit 1
fi

awk -v rtol="$rtol" -v atol="$atol" -v refname="$reference" -v outname="$output" '
function isnum(s) {
    return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eEdD][-+]?[0-9]+)?$/
}
function isnan(s) {
    return tolower(s) ~ /^[-+]?nan$/
}
function abs(x) { return x < 0 ? -x : x }
# One unit in the last printed digit of a numeric token (0 if it has no
# decimal point), so a rounding flip in the final printed digit is allowed.
function unit(s,    parts, n, mant, ex, nd) {
    n = split(tolower(s), parts, /[ed]/)
    mant = parts[1]
    ex = (n > 1) ? parts[2] + 0 : 0
    if (index(mant, ".") == 0) return 0
    nd = length(mant) - index(mant, ".")
    return 10 ^ (ex - nd)
}

# Pass 1: the reference file.
FNR == NR {
    nref = FNR
    ref_nf[FNR] = NF
    for (i = 1; i <= NF; ++i) {
        ref[FNR, i] = $i
        if (isnum($i)) {
            v = abs($i + 0)
            if (v > scale[i]) scale[i] = v
            if (v > filescale) filescale = v
        }
    }
    next
}

# Pass 2: the output file.
{
    nout = FNR
    if (FNR > nref) next
    if (NF != ref_nf[FNR]) {
        printf "%s:%d: %d fields, reference %s has %d\n", outname, FNR, NF, refname, ref_nf[FNR]
        failures++
        next
    }
    for (i = 1; i <= NF; ++i) {
        a = ref[FNR, i]
        b = $i
        if (isnum(a) && isnum(b)) {
            x = a + 0; y = b + 0
            diff = abs(x - y)
            m = abs(x) > abs(y) ? abs(x) : abs(y)
            s = scale[i] > 0 ? scale[i] : filescale
            ua = unit(a); ub = unit(b)
            tol = rtol * m + atol * s + (ua > ub ? ua : ub)
            if (s > 0 && diff / s > maxrel) { maxrel = diff / s; maxrel_at = FNR ":" i }
            if (diff > tol) {
                if (failures < 20)
                    printf "%s:%d field %d: %s differs from reference %s (|diff| = %g > %g)\n", \
                        outname, FNR, i, b, a, diff, tol
                failures++
            }
        } else if (isnan(a) && isnan(b)) {
            continue
        } else if (a != b) {
            if (failures < 20)
                printf "%s:%d field %d: \"%s\" differs from reference \"%s\"\n", outname, FNR, i, b, a
            failures++
        }
    }
}

END {
    if (nout != nref) {
        printf "%s: %d lines, reference %s has %d\n", outname, nout, refname, nref
        failures++
    }
    printf "%s: max deviation relative to column scale %g", outname, maxrel + 0
    if (maxrel_at != "") printf " (line:field %s)", maxrel_at
    printf "\n"
    if (failures) {
        printf "%s: %d mismatches against %s\n", outname, failures, refname
        exit 1
    }
}
' "$reference" "$output"
