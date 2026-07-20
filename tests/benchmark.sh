#!/usr/bin/env bash
# Benchmark script for ihsnpeaks parallel sweep validation.
# Usage: ./tests/benchmark.sh [--save-ref] [fmax] [degrees]
#   --save-ref  save current j=1 peaks as reference files, then exit
#   fmax        upper frequency bound (default: 100)
#   degrees     space-separated list of degrees (default: "1 8")
#
# Requires: bash, awk, the ihsnpeaks binary in project root.
# Target: OGLE-BLAP-035 (2240 measurements, ~3500d timebase).
# Reference files: tests/data/ref_<grid>_d<degree>.peaks (peaks only, no spectrum).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BINARY="$PROJECT_DIR/ihsnpeaks"
TARGET="$SCRIPT_DIR/data/OGLE-BLAP-035.dat"

FMAX="${1:-100}"
DEGREES="${2:-1 8}"
JOBS_LIST="1 2 4 8 16 24 32"
SAVE_REF=""

if [[ "${1:-}" == "--save-ref" ]]; then
    SAVE_REF="--save-ref"
    FMAX="${2:-100}"
    DEGREES="${3:-1 8}"
fi

if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: binary not found at $BINARY (run 'make' first)" >&2
    exit 1
fi

if [[ ! -f "$TARGET" ]]; then
    echo "ERROR: test data not found at $TARGET" >&2
    exit 1
fi

# time_run <jobs> <grid> <degree>
# Prints elapsed time in milliseconds (wall clock, peaks-only mode, no spectrum output).
time_run() {
    local jobs=$1 grid=$2 degree=$3
    local start end
    start=$(date +%s%N)
    "$BINARY" "$TARGET" "$FMAX" -g "$grid" -d "$degree" -j "$jobs" > /dev/null 2>&1
    end=$(date +%s%N)
    echo $(( (end - start) / 1000000 ))
}

# extract_peaks <grid> <degree>
# Runs j=1 and extracts sorted peaks (freq R²) — used for reference generation and comparison.
extract_peaks() {
    local grid=$1 degree=$2
    "$BINARY" "$TARGET" "$FMAX" -g "$grid" -d "$degree" -j 1 2>/dev/null | LC_ALL=C awk 'NR>3 && NF>=4 {printf "%.6f %.4f\n", $1, $4}'
}

# correctness_check <grid> <degree>
# Compares current j=1 output against saved reference (tests/data/ref_<grid>_d<degree>.peaks).
# If --save-ref is passed as $3, writes the reference file instead.
correctness_check() {
    local grid=$1 degree=$2 save_ref="${3:-}"
    local reffile="$SCRIPT_DIR/data/ref_${grid}_d${degree}.peaks"
    local current
    current=$(extract_peaks "$grid" "$degree")

    if [[ "$save_ref" == "--save-ref" ]]; then
        echo "$current" > "$reffile"
        echo "SAVED ($reffile)"
        return
    fi

    if [[ ! -f "$reffile" ]]; then
        echo "NO REF (run with --save-ref to generate)"
        return
    fi

    local expected
    expected=$(cat "$reffile")
    if [[ "$current" == "$expected" ]]; then
        echo "PASS"
    else
        echo "FAIL"
        diff <(echo "$expected") <(echo "$current") | head -10 | sed 's/^/  /'
    fi
}

echo "================================================================"
echo " ihsnpeaks benchmark"
echo " target: OGLE-BLAP-035 (2240 pts, ~3500d timebase)"
echo " fmax: $FMAX"
echo "================================================================"

for degree in $DEGREES; do
    for grid in ihs aov; do
        echo ""
        echo "--- $grid d=$degree ---"

        # Correctness
        result=$(correctness_check "$grid" "$degree" "$SAVE_REF")
        echo "  correctness: $result"

        # Skip benchmarks when saving references
        if [[ -n "$SAVE_REF" ]]; then
            continue
        fi

        # Benchmark
        baseline=""
        printf "  %4s  %8s  %8s\n" "j" "ms" "speedup"
        for jobs in $JOBS_LIST; do
            ms=$(time_run "$jobs" "$grid" "$degree")
            if [[ -z "$baseline" ]]; then
                baseline=$ms
            fi
            if [[ "$baseline" -gt 0 ]]; then
                speedup=$(LC_ALL=C awk "BEGIN {printf \"%.2f\", $baseline / $ms}")
            else
                speedup="N/A"
            fi
            printf "  %4d  %6d ms  %6sx\n" "$jobs" "$ms" "$speedup"
        done
    done
done

echo ""
echo "Done."
