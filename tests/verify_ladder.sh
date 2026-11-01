#!/usr/bin/env bash
# Compile and run the twiddle-ladder / vector-trig regression tests.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CC="${CC:-cc}"
OUT="${TMPDIR:-/tmp}/ihsnpeaks_ladder_test"

"$CC" -O2 -std=gnu23 -D_GNU_SOURCE \
    -I"$PROJECT_DIR/src" \
    -I"$PROJECT_DIR/src/nufft" \
    -I"$PROJECT_DIR/build/native/system" \
    -I"$PROJECT_DIR/include" \
    "$SCRIPT_DIR/ladder_test.c" \
    "$PROJECT_DIR/src/nufft/nufft1.c" \
    -o "$OUT" -lm

"$OUT"
