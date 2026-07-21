#!/bin/bash
# Verify parallel IHS/AoV determinism by running N times and comparing output.
# Usage: bash tests/verify_determinism.sh [binary] [data] [fmax] [runs] [extra_args...]
# Example:
#   bash tests/verify_determinism.sh ./ihsnpeaks tests/data/OGLE-BLAP-035.dat 50 10
#   bash tests/verify_determinism.sh ./ihsnpeaks tests/data/OGLE-BLAP-035.dat 50 10 -g aov

BINARY="${1:-./ihsnpeaks}"
DATA="${2:-tests/data/OGLE-BLAP-035.dat}"
FMAX="${3:-50}"
RUNS="${4:-10}"
shift 4 2>/dev/null
EXTRA_ARGS="$@"

echo "Binary: $BINARY"
echo "Data:   $DATA"
echo "Fmax:   $FMAX"
echo "Runs:   $RUNS"
echo "Extra:  $EXTRA_ARGS"
echo "---"

REFERENCE=$($BINARY "$DATA" "$FMAX" $EXTRA_ARGS 2>/dev/null)
if [ -z "$REFERENCE" ]; then
    echo "ERROR: No output from first run"
    exit 1
fi

PASS=true
for i in $(seq 2 $RUNS); do
    OUTPUT=$($BINARY "$DATA" "$FMAX" $EXTRA_ARGS 2>/dev/null)
    if [ "$OUTPUT" != "$REFERENCE" ]; then
        echo "FAIL: Run $i differs from run 1"
        diff <(echo "$REFERENCE") <(echo "$OUTPUT") | head -20
        PASS=false
    fi
done

if $PASS; then
    echo "PASS: All $RUNS runs produced identical output"
else
    echo "FAILED: Non-deterministic output detected"
    exit 1
fi
