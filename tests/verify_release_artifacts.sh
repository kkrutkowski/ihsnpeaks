#!/usr/bin/env bash
# Smoke-test the `make release-full` artifacts: the static x86-64 and ARM64
# Linux binaries and the universal macOS binary.
#
# x86-64 runs natively.  The ARM64 binary runs under qemu-aarch64 and the
# macOS binary under Darling when those emulators are installed; otherwise that
# check is skipped with a note, so the script still works on a bare host.  Every
# check is one 2240-point light curve on a 10 cycles/day grid with a single
# thread, which keeps the emulated runs to well under a second each.
#
# Usage: tests/verify_release_artifacts.sh [--require-emulators]
#   --require-emulators  treat a missing or unusable qemu-aarch64 / Darling as a
#                        failure instead of a skip (use this on CI images that
#                        ship both emulators)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA="tests/data/OGLE-BLAP-035.dat"
RUN_ARGS=(10 -j 1 -b 0 -n 3 -t 0)

REQUIRE_EMULATORS=0
for arg in "$@"; do
    case "$arg" in
        --require-emulators) REQUIRE_EMULATORS=1 ;;
        *)
            echo "unknown argument: $arg" >&2
            exit 2
            ;;
    esac
done

# Darling only sees the project tree and only through a path relative to the
# current directory, so every runner below works from the project root.
cd "$PROJECT_DIR"

DIST_DIR="${DIST_DIR:-dist}"
if [[ "$DIST_DIR" = /* ]]; then
    DIST_DIR="$(realpath --relative-to="$PROJECT_DIR" "$DIST_DIR")"
fi
if [[ "$DIST_DIR" = ../* ]]; then
    echo "DIST_DIR must live inside $PROJECT_DIR (Darling cannot reach absolute host paths)" >&2
    exit 2
fi

X86_BIN="$DIST_DIR/ihsnpeaks-linux-x86-64"
ARM_BIN="$DIST_DIR/ihsnpeaks-linux-arm64"
MACOS_BIN="$DIST_DIR/ihsnpeaks-macos"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ihsnpeaks_release.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

checks=0
skips=0
failed=0

pass() {
    checks=$((checks + 1))
    echo "ok   $*"
}

fail() {
    checks=$((checks + 1))
    failed=$((failed + 1))
    echo "FAIL $*" >&2
}

skip() {
    skips=$((skips + 1))
    echo "skip $*"
}

# show FILE: indent a tool's diagnostics under the current check.
show() {
    head -n 4 "$1" | sed 's/^/     /' >&2 || true
}

find_tool() {
    local candidate
    for candidate in "$@"; do
        if command -v "$candidate" > /dev/null 2>&1; then
            command -v "$candidate"
            return 0
        fi
    done
    return 1
}

# run_capture LABEL OUT CMD...: run CMD with stdout+stderr captured in OUT and
# report a failure instead of aborting the script (set -e).
run_capture() {
    local label=$1 out=$2 status=0
    shift 2
    "$@" > "$out" 2>&1 || status=$?
    if ((status != 0)); then
        fail "$label exited with status $status"
        show "$out"
        return 1
    fi
    return 0
}

# --- artifacts ------------------------------------------------------------

for bin in "$X86_BIN" "$ARM_BIN" "$MACOS_BIN"; do
    if [[ -s "$bin" ]]; then
        pass "artifact exists: $bin"
    else
        fail "artifact missing or empty: $bin (run 'make release-full' first)"
    fi
done
if ((failed)); then
    echo "verify_release_artifacts: FAIL (missing artifacts)" >&2
    exit 1
fi

X86_FILE="$(file -b "$X86_BIN")"
ARM_FILE="$(file -b "$ARM_BIN")"
MACOS_FILE="$(file -b "$MACOS_BIN")"

if [[ "$X86_FILE" == *"x86-64"* && "$X86_FILE" == *ELF* && "$X86_FILE" == *static* ]]; then
    pass "x86-64 artifact is a statically linked ELF"
else
    fail "x86-64 artifact is not a static x86-64 ELF: $X86_FILE"
fi
if [[ "$ARM_FILE" == *"aarch64"* && "$ARM_FILE" == *ELF* && "$ARM_FILE" == *static* ]]; then
    pass "ARM64 artifact is a statically linked ELF"
else
    fail "ARM64 artifact is not a static aarch64 ELF: $ARM_FILE"
fi
if [[ "$MACOS_FILE" == *"Mach-O universal"* && "$MACOS_FILE" == *x86_64* && "$MACOS_FILE" == *arm64* ]]; then
    pass "macOS artifact is a universal Mach-O (x86_64 + arm64)"
else
    fail "macOS artifact is not a universal Mach-O: $MACOS_FILE"
fi

EXPECTED_VERSION="$(sed -n 's/^#define IHSNPEAKS_VERSION "\(.*\)"/\1/p' src/main.c)"
if [[ -z "$EXPECTED_VERSION" ]]; then
    echo "cannot read IHSNPEAKS_VERSION from src/main.c" >&2
    exit 1
fi

# --- helpers --------------------------------------------------------------

# check_version LABEL OUT: the first line of --help must name the version that
# src/main.c declares.
check_version() {
    local label=$1 out=$2 reported
    reported="$(head -n 1 "$out")"
    if [[ "$reported" == "ihsnpeaks $EXPECTED_VERSION" ]]; then
        pass "$label reports $EXPECTED_VERSION"
    else
        fail "$label reports '$reported', expected 'ihsnpeaks $EXPECTED_VERSION'"
    fi
}

# extract_peaks RAW OUT: keep only the peak rows (leading whitespace + digit).
extract_peaks() {
    grep -E '^[[:space:]]+[0-9]' "$1" > "$2" || true
}

# compare_peaks REF CAND LABEL: CAND must emit the same number of peaks with the
# same frequencies, NLL values and amplitudes.  Different ISA levels accumulate
# in a different order, so the numeric comparison is tolerant.
compare_peaks() {
    local ref=$1 cand=$2 label=$3
    local ref_lines cand_lines
    ref_lines=$(wc -l < "$ref")
    cand_lines=$(wc -l < "$cand")
    if [[ "$ref_lines" -eq 0 ]]; then
        fail "$label produced no peaks and the x86-64 reference produced none either"
        return
    fi
    if [[ "$ref_lines" -ne "$cand_lines" ]]; then
        fail "$label produced $cand_lines peaks, x86-64 reference produced $ref_lines"
        show "$cand"
        return
    fi
    if awk '
        FNR == NR { freq[FNR] = $1; nll[FNR] = $2; amp[FNR] = $3; next }
        {
            fd = $1 - freq[FNR]; if (fd < 0) fd = -fd
            nd = $2 - nll[FNR]; if (nd < 0) nd = -nd
            ad = $3 - amp[FNR]; if (ad < 0) ad = -ad
            if (fd > 5e-5 || nd > 0.02 || ad > 2e-3) {
                printf "     row %d: got (%.6f, %.2f, %.3f), expected (%.6f, %.2f, %.3f)\n",
                    FNR, $1, $2, $3, freq[FNR], nll[FNR], amp[FNR]
                bad = 1
            }
        }
        END { exit bad }
    ' "$ref" "$cand"; then
        pass "$label matches the x86-64 peaks"
    else
        fail "$label peaks differ from the x86-64 reference"
    fi
}

# preflight_runner LABEL CMD...: returns 0 when the runner can start, 1 when it
# is unavailable, 2 when it is installed but broken.
preflight_runner() {
    local label=$1 status=0
    shift
    "$@" > /dev/null 2>&1 || status=$?
    if ((status == 0)); then
        return 0
    fi
    if ((REQUIRE_EMULATORS)); then
        echo "     $label is installed but could not start (command exited with status $status)" >&2
        return 2
    fi
    return 1
}

emulator_verdict() {
    local label=$1 status=$2
    if ((status == 0)); then
        return 0
    fi
    if ((status == 2)); then
        fail "$label is installed but cannot run the artifacts"
        return 1
    fi
    skip "$label not available on this host, skipping its checks"
    return 1
}

# --- x86-64 (reference) ---------------------------------------------------

if run_capture "x86-64 --help" "$WORK/x86.help" "$X86_BIN" --help; then
    check_version "x86-64" "$WORK/x86.help"
fi

if run_capture "x86-64 dispatch list" "$WORK/x86.dispatch" env IHSNPEAKS_DISPATCH=list "$X86_BIN"; then
    missing=()
    for variant in x86-64 x86-64-v2 x86-64-v2+avx x86-64-v3 x86-64-v4; do
        grep -qxF "  $variant" "$WORK/x86.dispatch" || missing+=("$variant")
    done
    if ((${#missing[@]} == 0)); then
        pass "x86-64 artifact dispatches all 5 variants"
    else
        fail "x86-64 artifact is missing dispatch variant(s): ${missing[*]}"
    fi
fi

: > "$WORK/x86.peaks"
if run_capture "x86-64 periodogram" "$WORK/x86.run" "$X86_BIN" "$DATA" "${RUN_ARGS[@]}"; then
    extract_peaks "$WORK/x86.run" "$WORK/x86.peaks"
    if [[ -s "$WORK/x86.peaks" ]]; then
        pass "x86-64 artifact computes a periodogram ($(wc -l < "$WORK/x86.peaks") peaks)"
    else
        fail "x86-64 artifact produced no peaks on $DATA"
    fi
fi

# --- ARM64 under qemu -----------------------------------------------------

QEMU="${QEMU_AARCH64:-$(find_tool qemu-aarch64 qemu-aarch64-static || true)}"
if [[ -z "$QEMU" ]]; then
    if ((REQUIRE_EMULATORS)); then
        fail "qemu-aarch64 not found, cannot test the ARM64 artifact"
    else
        skip "qemu-aarch64 not installed, skipping the ARM64 checks"
    fi
else
    arm_status=0
    preflight_runner "$QEMU" "$QEMU" --version || arm_status=$?
    if emulator_verdict "$QEMU" "$arm_status"; then
        if run_capture "ARM64 --help" "$WORK/arm.help" "$QEMU" "$ARM_BIN" --help; then
            check_version "ARM64 (qemu)" "$WORK/arm.help"
        fi

        if run_capture "ARM64 dispatch list" "$WORK/arm.dispatch" env IHSNPEAKS_DISPATCH=list "$QEMU" "$ARM_BIN"; then
            if grep -qxF '  arm64-neon' "$WORK/arm.dispatch" && grep -qE '^  arm64-sve' "$WORK/arm.dispatch"; then
                pass "ARM64 artifact dispatches the NEON and SVE variants"
            else
                fail "ARM64 artifact does not advertise the expected NEON/SVE variants"
                show "$WORK/arm.dispatch"
            fi
        fi

        if run_capture "ARM64 periodogram" "$WORK/arm.run" "$QEMU" "$ARM_BIN" "$DATA" "${RUN_ARGS[@]}"; then
            extract_peaks "$WORK/arm.run" "$WORK/arm.peaks"
            compare_peaks "$WORK/x86.peaks" "$WORK/arm.peaks" "ARM64 (qemu)"
        fi
    fi
fi

# --- macOS under Darling --------------------------------------------------

DARLING="${DARLING:-$(find_tool darling || true)}"
if [[ -z "$DARLING" ]]; then
    if ((REQUIRE_EMULATORS)); then
        fail "Darling not found, cannot test the macOS artifact"
    else
        skip "Darling not installed, skipping the macOS checks"
    fi
else
    macos_status=0
    preflight_runner "Darling" "$DARLING" shell /usr/bin/true || macos_status=$?
    if emulator_verdict "Darling" "$macos_status"; then
        if run_capture "macOS --help" "$WORK/macos.help" "$DARLING" shell "$MACOS_BIN" --help; then
            check_version "macOS (Darling)" "$WORK/macos.help"
        fi

        if run_capture "macOS periodogram" "$WORK/macos.run" "$DARLING" shell "$MACOS_BIN" "$DATA" "${RUN_ARGS[@]}"; then
            extract_peaks "$WORK/macos.run" "$WORK/macos.peaks"
            compare_peaks "$WORK/x86.peaks" "$WORK/macos.peaks" "macOS (Darling)"
        fi
    fi
fi

# --- summary --------------------------------------------------------------

if ((failed)); then
    echo "verify_release_artifacts: FAIL ($failed of $checks checks failed, $skips skipped)" >&2
    exit 1
fi
echo "verify_release_artifacts: PASS ($checks checks, $skips skipped)"
