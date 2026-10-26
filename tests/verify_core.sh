#!/usr/bin/env bash
# Compile and run the core C regression tests for the d=1 condition/Cybenko
# path, the forward-only Levinson solver, and the vectorized Wilkie correction.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CC="${CC:-cc}"
OUT="${TMPDIR:-/tmp}/ihsnpeaks_core_fixes_test"

"$CC" -O2 -std=gnu23 -D_GNU_SOURCE \
    -I"$PROJECT_DIR/src" \
    -I"$PROJECT_DIR/src/nufft" \
    -I"$PROJECT_DIR/build/native/system" \
    -I"$PROJECT_DIR/include" \
    "$SCRIPT_DIR/core_fixes_test.c" \
    -o "$OUT" -lm

"$OUT"
