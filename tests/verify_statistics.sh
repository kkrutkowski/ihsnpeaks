#!/usr/bin/env bash
# CLI-level regression checks for the statistic defaults, IHS corrected mode,
# validation errors, and corrected parallel determinism.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BINARY="${1:-$PROJECT_DIR/ihsnpeaks}"
DATA="${2:-$PROJECT_DIR/tests/data/OGLE-BLAP-035.dat}"

if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: binary not found at $BINARY (run 'make' first)" >&2
    exit 1
fi
if [[ ! -f "$DATA" ]]; then
    echo "ERROR: data file not found at $DATA" >&2
    exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$DATA" "$TMP/lc.dat"

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

check_default() {
    local d=$1 method=$2 stat=$3
    local out
    out="$("$BINARY" "$TMP/lc.dat" 50 -d "$d" --debug 2>&1)"
    grep -q "Periodogram method: $method" <<<"$out" || fail "d=$d method != $method"
    grep -q "Statistic: $stat" <<<"$out" || fail "d=$d statistic != $stat"
    echo "PASS default d=$d -> $method/$stat"
}

check_default 1 ihs raw
check_default 4 ihs raw
check_default 5 aov raw
check_default 7 aov raw
check_default 8 aov bayes
check_default 9 aov bayes

out="$("$BINARY" "$TMP/lc.dat" 50 -g ihs -d 9 --debug 2>&1)"
grep -q "Periodogram method: ihs" <<<"$out" || fail "explicit ihs d=9 method"
grep -q "Statistic: raw" <<<"$out" || fail "explicit ihs d=9 statistic"

expect_fail() {
    local label=$1
    shift
    if "$BINARY" "$TMP/lc.dat" 50 "$@" >/dev/null 2>&1; then
        fail "$label should have failed"
    fi
    echo "PASS rejected: $label"
}

expect_fail "corrected with AoV" -g aov -d 5 --statistic corrected
expect_fail "bayes with IHS" -g ihs -d 3 --statistic bayes
expect_fail "corrected with direct-eval mode" -g ihs -d 3 --statistic corrected -m 5

"$BINARY" "$TMP/lc.dat" 50 -g ihs -d 3 --statistic raw --ddet 0 --epsilon 0 -s -n 0 -j 1 >/dev/null 2>&1
mv "$TMP/lc.tsv" "$TMP/raw.tsv"
"$BINARY" "$TMP/lc.dat" 50 -g ihs -d 3 --statistic corrected --ddet 0 --epsilon 0 -s -n 0 -j 1 >/dev/null 2>&1
mv "$TMP/lc.tsv" "$TMP/corrected.tsv"

python3 - "$TMP/raw.tsv" "$TMP/corrected.tsv" <<'PY'
import sys
import math

raw = []
corr = []
for path, dest in ((sys.argv[1], raw), (sys.argv[2], corr)):
    with open(path, "r", encoding="utf-8") as src:
        for line in src:
            parts = line.split()
            if len(parts) < 2:
                continue
            value = float(parts[1])
            if not math.isfinite(value) or value < 0.0:
                raise SystemExit(f"non-finite/negative spectrum value in {path}: {value}")
            dest.append(value)

if len(raw) != len(corr) or not raw:
    raise SystemExit("raw/corrected spectrum length mismatch or empty output")

max_diff = max(abs(a - b) for a, b in zip(raw, corr))
if max_diff <= 1.0e-6:
    raise SystemExit("corrected spectrum is identical to raw")
print(f"PASS corrected spectrum differs from raw (max |diff| = {max_diff:.3g})")
PY

"$BINARY" "$TMP/lc.dat" 50 -g ihs -d 3 --statistic corrected -p -n 3 -j 1 >/dev/null 2>&1
echo "PASS corrected prewhiten run"

# Determinism is checked per thread count.  Serial vs sliced-parallel runs
# differ by small floating-point amounts by design (see process.h), so this
# mirrors tests/verify_determinism.sh rather than requiring -j 1 == -j 8.
for j in 1 8; do
    "$BINARY" "$TMP/lc.dat" 50 -g ihs -d 3 --statistic corrected -n 10 -c 10 -j "$j" >"$TMP/run_a.txt" 2>&1
    "$BINARY" "$TMP/lc.dat" 50 -g ihs -d 3 --statistic corrected -n 10 -c 10 -j "$j" >"$TMP/run_b.txt" 2>&1
    cmp -s "$TMP/run_a.txt" "$TMP/run_b.txt" || fail "corrected output is not deterministic at -j $j"
done
echo "PASS corrected determinism at -j 1 and -j 8"

echo "verify_statistics.sh: PASS"
