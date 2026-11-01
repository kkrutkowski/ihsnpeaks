#!/usr/bin/env bash
# Functional tests for .tsv list input and --prefix.
#
# Usage: ./tests/verify_list_input.sh [path-to-binary]
# Requires: bash, awk/sed/cut, bash builtins; the test photometry in tests/data.
#
# List mode is always batch: every run below writes <list>.out.tsv next to the
# list argument and never modifies the list itself. Runs use -b 0 -j 1 so row
# order is deterministic (lists are processed in file order).

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BINARY="${1:-$PROJECT_DIR/ihsnpeaks}"
SOURCE_DAT="$SCRIPT_DIR/data/OGLE-BLAP-035.dat"

if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: binary not found or not executable: $BINARY (run 'make' first)" >&2
    exit 1
fi
if [[ ! -f "$SOURCE_DAT" ]]; then
    echo "ERROR: test photometry not found: $SOURCE_DAT" >&2
    exit 1
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ihsnpeaks_list.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

PASS=0
FAIL=0
FAILURES=()

pass() {
    PASS=$((PASS + 1))
    echo "  ok: $1"
}

fail() {
    FAIL=$((FAIL + 1))
    FAILURES+=("$1")
    echo "  FAIL: $1" >&2
}

LAST_OUT="$WORK/last_stdout.txt"
LAST_ERR="$WORK/last_stderr.txt"
LAST_EXIT=0

run() {
    "$BINARY" "$@" >"$LAST_OUT" 2>"$LAST_ERR"
    LAST_EXIT=$?
}

expect_exit() {
    if [[ "$LAST_EXIT" -eq "$1" ]]; then
        pass "$2"
    else
        fail "$2 (expected exit $1, got $LAST_EXIT; stderr: $(head -3 "$LAST_ERR"))"
    fi
}

expect_file() {
    if [[ -f "$1" ]]; then
        pass "$2"
    else
        fail "$2 (missing $1)"
    fi
}

expect_no_file() {
    if [[ -e "$1" ]]; then
        fail "$2 (unexpected $1)"
    else
        pass "$2"
    fi
}

# expect_labels <output.tsv> <label...>
expect_labels() {
    local file="$1"
    shift
    local expected got
    expected="$(printf '%s\n' "$@")"
    got="$(tail -n +2 "$file" | cut -f1)"
    if [[ "$got" == "$expected" ]]; then
        pass "labels in $(basename "$file") match list order"
    else
        fail "labels in $file: got [$got], expected [$expected]"
    fi
}

expect_stderr_contains() {
    local needle="$1" name="$2"
    if grep -qF -- "$needle" "$LAST_ERR"; then
        pass "$name"
    else
        fail "$name (stderr did not contain '$needle': $(head -3 "$LAST_ERR"))"
    fi
}

expect_same_bytes() {
    if cmp -s "$1" "$2"; then
        pass "$3"
    else
        fail "$3 ($1 differs from $2)"
    fi
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
mkdir -p "$WORK/data" "$WORK/elsewhere" "$WORK/lists" "$WORK/onefile"
cp "$SOURCE_DAT" "$WORK/data/a.dat"
cp "$SOURCE_DAT" "$WORK/data/b.dat"
cp "$SOURCE_DAT" "$WORK/data/c.dat"
cp "$SOURCE_DAT" "$WORK/onefile/only.dat"

cat >"$WORK/list.tsv" <<'EOF'
Input_file	unused	columns
a.dat	0.1	0.2
# comment line, not a file

b.dat	nonsense
EOF
cp "$WORK/list.tsv" "$WORK/list.tsv.bak"

FLAGS=(-m 0 -t 0 -b 0 -j 1)
cd "$WORK"

echo "ihsnpeaks list-input tests"
echo "  binary: $BINARY"
echo "  work:   $WORK"

# ---------------------------------------------------------------------------
# 1. Basic list: header, extra columns, comment and blank line are ignored
# ---------------------------------------------------------------------------
echo "1. basic list parsing and output naming"
run "$WORK/list.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "list run exits 0"
expect_file "$WORK/list.out.tsv" "results written to <list>.out.tsv"
expect_same_bytes "$WORK/list.tsv" "$WORK/list.tsv.bak" "input list is not modified"
expect_labels "$WORK/list.out.tsv" a.dat b.dat
if [[ "$(wc -l <"$WORK/list.out.tsv")" -eq 3 ]]; then
    pass "output has header + one row per listed file"
else
    fail "unexpected output row count: $(wc -l <"$WORK/list.out.tsv")"
fi

# ---------------------------------------------------------------------------
# 2. Prefix is applied to list entries, never to the list path itself
# ---------------------------------------------------------------------------
echo "2. prefix applies to entries, not to the list"
run "$WORK/lists/run.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 1 "missing list under a fresh name exits non-zero"
cp "$WORK/list.tsv" "$WORK/lists/run.tsv"
run "lists/run.tsv" 1.0 "${FLAGS[@]}" --prefix "data"
expect_exit 0 "relative list path + relative prefix exits 0"
expect_file "$WORK/lists/run.out.tsv" "output stays next to the literal list argument"
expect_no_file "$WORK/data/run.out.tsv" "output is not written into the prefixed data root"
expect_labels "$WORK/lists/run.out.tsv" a.dat b.dat

# ---------------------------------------------------------------------------
# 3. A previous output list can be fed back in (header skipped)
# ---------------------------------------------------------------------------
echo "3. output list is re-feedable"
cp "$WORK/list.out.tsv" "$WORK/chain.tsv"
run "$WORK/chain.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data/"
expect_exit 0 "re-feeding an .out.tsv exits 0"
expect_file "$WORK/chain.out.tsv" "chained output written"
expect_labels "$WORK/chain.out.tsv" a.dat b.dat

# ---------------------------------------------------------------------------
# 4. CRLF and a UTF-8 BOM are tolerated
# ---------------------------------------------------------------------------
echo "4. CRLF + BOM"
printf '\xEF\xBB\xBFa.dat\r\nb.dat\r\n' >"$WORK/crlf.tsv"
run "$WORK/crlf.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "CRLF/BOM list exits 0"
expect_labels "$WORK/crlf.out.tsv" a.dat b.dat

# ---------------------------------------------------------------------------
# 5. A list always writes its output file, even with one entry
# ---------------------------------------------------------------------------
echo "5. single-entry list"
printf 'a.dat\n' >"$WORK/single.tsv"
run "$WORK/single.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "single-entry list exits 0"
expect_file "$WORK/single.out.tsv" "single-entry list still writes <list>.out.tsv"
expect_labels "$WORK/single.out.tsv" a.dat

# ---------------------------------------------------------------------------
# 6. Absolute entries bypass the prefix (with a warning) and stay verbatim
# ---------------------------------------------------------------------------
echo "6. absolute list entries"
printf '%s\n' "$WORK/data/a.dat" >"$WORK/abs.tsv"
run "$WORK/abs.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "absolute entry exits 0"
expect_stderr_contains "--prefix ignored" "absolute entry reports the ignored prefix"
expect_labels "$WORK/abs.out.tsv" "$WORK/data/a.dat"

# ---------------------------------------------------------------------------
# 7. Invalid entries are skipped with warnings; all-invalid lists fail
# ---------------------------------------------------------------------------
echo "7. invalid entries"
printf 'a.dat\nmissing.dat\n' >"$WORK/missing.tsv"
run "$WORK/missing.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "list with one missing entry still exits 0"
expect_stderr_contains "skipping list entry 'missing.dat'" "missing entry is reported"
expect_stderr_contains "skipped 1 of 2" "skip summary is reported"
expect_labels "$WORK/missing.out.tsv" a.dat

echo "not a light curve" >"$WORK/data/notes.txt"
printf 'a.dat\nnotes.txt\n' >"$WORK/unsupported.tsv"
run "$WORK/unsupported.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "list with one unsupported entry still exits 0"
expect_stderr_contains "unsupported format" "unsupported entry is reported"
expect_labels "$WORK/unsupported.out.tsv" a.dat

: >"$WORK/data/empty.dat"
printf 'empty.dat\na.dat\n' >"$WORK/empty_phot.tsv"
run "$WORK/empty_phot.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "list with one empty photometry file still exits 0"
expect_stderr_contains "No data rows found" "empty photometry file is reported"
expect_labels "$WORK/empty_phot.out.tsv" a.dat

printf 'empty.dat\n' >"$WORK/all_empty_phot.tsv"
run "$WORK/all_empty_phot.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 1 "list of only empty photometry files exits non-zero"

: >"$WORK/empty.tsv"
run "$WORK/empty.tsv" 1.0 "${FLAGS[@]}"
expect_exit 1 "empty list exits non-zero"

printf 'Input_file\ta\tb\n' >"$WORK/header_only.tsv"
run "$WORK/header_only.tsv" 1.0 "${FLAGS[@]}"
expect_exit 1 "header-only list exits non-zero"

printf 'missing1.dat\nmissing2.dat\n' >"$WORK/all_bad.tsv"
run "$WORK/all_bad.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 1 "list with no usable entries exits non-zero"

run "$WORK/does-not-exist.tsv" 1.0 "${FLAGS[@]}"
expect_exit 1 "missing list exits non-zero"

# ---------------------------------------------------------------------------
# 8. -s writes spectra next to the resolved photometry, not next to the list
# ---------------------------------------------------------------------------
echo "8. -s with a prefixed list"
rm -f "$WORK/data/a.tsv" "$WORK/data/b.tsv"
run "$WORK/list.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data" -s
expect_exit 0 "list run with -s exits 0"
expect_file "$WORK/data/a.tsv" "spectrum written beside the prefixed photometry (a)"
expect_file "$WORK/data/b.tsv" "spectrum written beside the prefixed photometry (b)"
expect_same_bytes "$WORK/list.tsv" "$WORK/list.tsv.bak" "input list survives -s"
expect_labels "$WORK/list.out.tsv" a.dat b.dat

# ---------------------------------------------------------------------------
# 9. File and directory targets keep their old behavior; --prefix only
#    redirects where photometry is read from.
# ---------------------------------------------------------------------------
echo "9. file/directory targets with and without --prefix"
cp "$WORK/list.tsv" "$WORK/keep.tsv"
run "$WORK/keep.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_same_bytes "$WORK/keep.tsv" "$WORK/list.tsv.bak" "list copy still unmodified"

run "$WORK/data" 1.0 "${FLAGS[@]}"
expect_exit 0 "directory target exits 0"
expect_file "$WORK/data.tsv" "directory summary is written next to the literal target"
if [[ "$(wc -l <"$WORK/data.tsv")" -eq 4 ]]; then
    pass "directory summary has header + 3 rows"
else
    fail "unexpected directory row count: $(wc -l <"$WORK/data.tsv")"
fi

run "data" 1.0 "${FLAGS[@]}" --prefix "$WORK/"
expect_exit 0 "prefixed relative directory exits 0"
expect_file "$WORK/data.tsv" "prefixed directory summary follows the literal argument (cwd)"
cd "$WORK/elsewhere"
run "data" 1.0 "${FLAGS[@]}" --prefix "$WORK/"
cd "$WORK"
expect_exit 0 "prefixed directory from another cwd exits 0"
expect_file "$WORK/elsewhere/data.tsv" "prefixed directory summary stays in the invocation directory"
if cut -f1 "$WORK/elsewhere/data.tsv" | tail -n +2 | grep -q '^data/'; then
    pass "directory labels use the target as typed (prefix-free)"
else
    fail "directory labels unexpectedly contain the resolved prefix: $(cut -f1 "$WORK/elsewhere/data.tsv" | head -2)"
fi

# Single file: raw path does not exist locally, only prefix+target does.
cd "$WORK/elsewhere"
run "a.dat" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
cd "$WORK"
expect_exit 0 "single file resolved through --prefix exits 0"
if grep -q "^File: a.dat" "$LAST_OUT"; then
    pass "single-file label is the target as typed"
else
    fail "single-file terminal header did not use the prefix-free label"
fi

# A directory with a single file keeps its existing terminal-only behavior.
run "$WORK/onefile" 1.0 "${FLAGS[@]}"
expect_exit 0 "one-file directory exits 0"
expect_no_file "$WORK/onefile.tsv" "one-file directory still prints to the terminal only"

# ---------------------------------------------------------------------------
# 10. List mode must compute the same numbers as a directory batch.
# ---------------------------------------------------------------------------
echo "10. list vs directory equivalence"
printf 'a.dat\nb.dat\nc.dat\n' >"$WORK/equiv.tsv"
run "$WORK/equiv.tsv" 1.0 "${FLAGS[@]}" --prefix "$WORK/data"
expect_exit 0 "equivalence list run exits 0"
run "$WORK/data" 1.0 "${FLAGS[@]}"
expect_exit 0 "equivalence directory run exits 0"
if diff -q <(tail -n +2 "$WORK/equiv.out.tsv" | cut -f2- | sort) <(tail -n +2 "$WORK/data.tsv" | cut -f2- | sort) >/dev/null; then
    pass "list and directory runs produce identical peak columns"
else
    fail "list and directory peak columns differ"
fi

# ---------------------------------------------------------------------------
# 11. Help documents the new surface.
# ---------------------------------------------------------------------------
echo "11. help text"
run --help
expect_exit 0 "--help exits 0"
if grep -q -- "--prefix" "$LAST_OUT"; then
    pass "help documents --prefix"
else
    fail "help does not mention --prefix"
fi
if grep -qi "list" "$LAST_OUT"; then
    pass "help documents list input"
else
    fail "help does not mention list input"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
if [[ "$FAIL" -eq 0 ]]; then
    echo "PASS: $PASS checks"
    exit 0
fi
echo "FAIL: $FAIL of $((PASS + FAIL)) checks failed" >&2
for failure in "${FAILURES[@]}"; do
    echo "  - $failure" >&2
done
exit 1
