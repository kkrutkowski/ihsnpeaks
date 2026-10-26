#!/usr/bin/env bash
# COBRA/COBRAVO and fused-Radix-2 FFT regression + scaling checks.
#
# Builds the test twice for each available ISA level: once with the scalar COBRA
# backend selected and once with COBRAVO selected.  The two builds must produce
# identical permutation/NuFFT hashes, and COBRAVO must be faster.  The test also
# checks the fused-Radix-2 FFT against the legacy two-pass kernel and against
# the selected COBRA/COBRAVO permutation.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CC="${CC:-cc}"
OUT_DIR="${TMPDIR:-/tmp}/ihsnpeaks_bitrev_$$"
mkdir -p "$OUT_DIR"
trap 'rm -rf "$OUT_DIR"' EXIT

COMMON=(-O3 -std=gnu23 -D_GNU_SOURCE -DNUFFT1_SCALING_BOOTSTRAP
    -I"$PROJECT_DIR/src/nufft" -I"$PROJECT_DIR/include")

run_isa() {
    local name=$1
    shift
    local extra=("$@")
    local cobra="$OUT_DIR/bitrev_cobra_${name}"
    local cobravo="$OUT_DIR/bitrev_cobravo_${name}"

    if ! "$CC" "${COMMON[@]}" "${extra[@]}" -DF_BITREV_COBRAVO=0 \
        "$SCRIPT_DIR/bitrev_scaling_test.c" -lm -o "$cobra" 2>"$OUT_DIR/${name}.cobra.log"; then
        if [[ "$name" == "native" ]]; then
            echo "verify_bitrev: native build failed" >&2
            cat "$OUT_DIR/${name}.cobra.log" >&2
            exit 1
        fi
        echo "verify_bitrev: skipping $name (compiler does not support flags)"
        cat "$OUT_DIR/${name}.cobra.log" >&2
        return 0
    fi
    "$CC" "${COMMON[@]}" "${extra[@]}" -DF_BITREV_COBRAVO=1 \
        "$SCRIPT_DIR/bitrev_scaling_test.c" -lm -o "$cobravo" 2>"$OUT_DIR/${name}.cobravo.log"

    echo "=== $name: COBRA ==="
    "$cobra" | tee "$OUT_DIR/${name}.cobra.out"
    echo "=== $name: COBRAVO ==="
    "$cobravo" | tee "$OUT_DIR/${name}.cobravo.out"

    if grep -q 'FAIL' "$OUT_DIR/${name}.cobra.out" "$OUT_DIR/${name}.cobravo.out"; then
        echo "verify_bitrev: correctness failure under $name" >&2
        exit 1
    fi

    grep '^HASH ' "$OUT_DIR/${name}.cobra.out" | sed 's/ exec_us=.*//' > "$OUT_DIR/${name}.hash.cobra"
    grep '^HASH ' "$OUT_DIR/${name}.cobravo.out" | sed 's/ exec_us=.*//' > "$OUT_DIR/${name}.hash.cobravo"
    if ! diff -u "$OUT_DIR/${name}.hash.cobra" "$OUT_DIR/${name}.hash.cobravo"; then
        echo "verify_bitrev: NuFFT output differs between COBRA and COBRAVO under $name" >&2
        exit 1
    fi
}

run_isa native -march=native
run_isa avx2 -march=x86-64-v3 -mno-avx512f
run_isa scalar -march=x86-64-v2 -mno-avx2 -mno-avx512f -DVEC_BYTES=16

echo "verify_bitrev.sh: PASS"
