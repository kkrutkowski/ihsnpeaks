# ihsnpeaks — Project Guide

## Project Overview

`ihsnpeaks` is a high-performance C command-line periodogram utility for astronomical light-curve analysis, licensed under GPLv3. It is designed as a drop-in replacement for legacy tools (FNPEAKS/MPEAKS, aovdist, FastChi2).

**Key capabilities:**
- IHS/Rayleigh-style multi-harmonic periodograms
- AoV(MH), AoVMH(W), FastChi² periodogram implementations
- Gaussian blur-based Supersmoother-like periodogram (gbls/gbaw)
- BLS (Box Least Squares) transit search
- ANCOVA and F-test for inequality-of-variance reevaluation
- Batch processing of OGLE-format `.dat` photometric files
- Native multithreading (pthreads) with thread-pool for parallel batch processing
- Built-in PSWF-based NuFFT backend (no FFTW3 dependency)

**Architecture:**
- Single-file CLI binary: `src/main.c` is the entry point
- Core periodogram logic lives in `src/process.h` (large inline header)
- NuFFT backend in `src/nufft/nufft1.c` and `nufft1.h`
- SIMD vectorization via GCC/Clang vector extensions (`src/utils/simd.h`)
- Uses klib (kvec, kthread, ketopt, kopen) for data structures and helpers
- Uses sds (Redis dynamic strings) for string handling
- Runtime dispatch: single binary selects the optimal CPU feature variant at startup
- Static musl builds for Linux; universal Mach-O binary for macOS

**Repository structure:**

```
.
├── src/
│   ├── main.c           — Entry point, CLI argument parsing
│   ├── params.h         — Parameter struct, init, parsing, help text
│   ├── metadata.h       — File discovery, NuFFT plan sizing, plan cache init
│   ├── process.h        — Core periodogram evaluation pipeline (780+ lines, main logic)
│   ├── profile.h        — Instrumentation/profiling macros (IHSNPEAKS_PROFILE)
│   ├── nufft/           — PSWF-based non-uniform FFT backend
│   │   ├── nufft1.c     — NuFFT implementation
│   │   ├── nufft1.h     — NuFFT public API
│   │   ├── scaling.c    — PSWF scaling coefficients generator
│   ├── utils/           — Shared utilities
│   │   ├── common.h     — Core types: buffer_t, parameters, peak_t, enums
│   │   ├── simd.h       — GCC vector extension SIMD wrappers (ln_ps, correctPower)
│   │   ├── aov.h        — AoV-based periodogram helpers
│   │   ├── convolution.h — Gaussian blur / BLS direct evaluation
│   │   ├── readout.h    — Buffer allocation, I/O, peak output formatting
│   │   ├── compat.h     — C11 aligned_alloc fallback for C99
│   │   ├── trig.h       — Trig utilities (sin2pif_tls, cos2pif_tls)
│   │   ├── spectrum.h   — Spectrum column capture
│   │   └── strings.h    — String formatting helpers
│   └── debug/           — Debugging helpers (not part of release)
│       ├── photview     — Python photometry viewer
│       └── spec_viewer.py — Python spectrum viewer
├── include/
│   ├── klib/            — External klib headers (kvec, kthread, ketopt, kopen)
│   ├── qfits/           — qfits library header
│   ├── fast_convert.h   — Fast float parsing
│   ├── fdist.h          — F-distribution CDF
│   └── sds.h            — Redis SDS dynamic strings
├── dispatch/
│   ├── variants.py      — Runtime dispatch variant definitions
│   ├── build_release.py — x86-64 musl release builder
│   ├── build_release_arm64.py — ARM64 musl release builder
│   ├── build_release_macos.py — macOS universal binary builder
│   └── zig_musl_cc.sh   — Zig-based cross-compilation wrapper
├── Dockerfile.release   — Multi-stage Docker build for release artifacts
├── makefile             — Primary build system (GNU Make)
└── README.md            — User-facing documentation
```

## Building

### Native (local machine)

```sh
make
# or equivalently:
make native
```

Uses `-march=native` and the host compiler (gcc, clang, or icx). Outputs the `ihsnpeaks` binary.

- Auto-detects C standard: prefers `-std=gnu23`, falls back to `gnu11`, then `gnu99`.
- Supports optional mimalloc allocator: `make MIMALLOC=1`.
- Runs `strip -s` on the final binary.

### Release with runtime dispatch (Docker)

```sh
make release-linux      # x86-64 + ARM64 Linux musl binaries
make release-full       # Linux + macOS universal binary
make release-x86        # x86-64 Linux only
make release-arm        # ARM64 Linux only
make release-macos      # macOS universal binary only
```

Linux artifacts use Docker (Alpine 3.23) with zig cross-compilation for musl static linking.
The macOS build uses `zig cc` + `llvm-lipo` to create a universal Mach-O binary.

### Profile build

```sh
make profile         # Builds a profiled binary with -IHSNPEAKS_PROFILE=1
make profile-run     # Runs profiling across thread counts (1, 2, 4, 8, 16, 32)
make profile-clean   # Removes profile artifacts
```

Configure via `PROFILE_ISA` (e.g., `make profile PROFILE_ISA=x86-64-v3`), `PROFILE_TARGET`, `PROFILE_FMAX`, and `PROFILE_JOBS`.

### Other make targets

```sh
make install              # cp ihsnpeaks /usr/local/bin/ (requires root)
make clean                # rm -rf build/
make format               # clang-format all .c/.h files in src/ and include/
make check_compiler       # Print detected compiler, C standard, CPU features
```

## Running

```
ihsnpeaks <target> <fmax> [options]
```

- `target`: path to a `.dat` file or a directory of `.dat` files
- `fmax`: upper frequency bound for the periodogram grid

**Common options:**

| Flag | Description |
|------|-------------|
| `-d N` | Number of harmonics/terms (default: 3) |
| `-m N` | Peak evaluation mode 0–6 (default: 2) |
| `-g METHOD` | Periodogram method: `ihs` (default), `aov`, `aovmh`, `aobmhw`, `chi`, `chi2`, `fastchi2` |
| `-e METHOD` | Peak evaluation: `gbls[alpha]`, `gbaw[alpha]`, or `bls[a,b,count]` |
| `-t THRESHOLD` | Detection threshold in dB (default: 10.0) |
| `-n N` | Max number of peaks to report (default: 10) |
| `-s` | Save full spectrum to `.tsv` files |
| `-p` | Prewhiten (subtract detected signals before finding next peaks) |
| `-j N` | Limit worker threads (default: 0 = all available) |
| `-o N` | Oversampling factor (default: 5.0) |
| `-f FMIN` | Lower frequency bound (default: 2/delta_t) |
| `--period` | Output periods instead of frequencies |
| `--nufft MODE` | NuFFT backend: `43`/`pswf43` (default) or `21`/`pswf21` |
| `--debug` | Print parsed parameters before computation |
| `-h` | Print help |

**Useful environment variables:**
- `IHSNPEAKS_DISPATCH=<variant>` — force a specific dispatch variant for testing
- `CC` — select a specific C compiler (default: `cc`)

## Development Conventions

- **Language:** GNU C (auto-detected: C23 > C11 > C99), `_GNU_SOURCE` always defined.
- **Style:** `clang-format` with `.clang-format` config at repo root and in subdirectories. Run `make format` before committing.
- **Headers as implementation units:** Core logic uses static inline functions in `.h` files (process.h, params.h, metadata.h) — no separate `.c` files for these modules. The NuFFT backend is the main exception (nufft1.c + nufft1.h).
- **No dynamic dispatch dispatcher source:** Runtime dispatch resolves via `IHSNPEAKS_DISPATCH` or is baked into release builds at link time.
- **Alignment:** All buffers use 64-byte aligned allocation (`aligned_alloc(64, ...)`). A C99 fallback keeps aligned allocations even on pre-C11 systems.
- **SIMD abstraction:** Vector code uses GCC/Clang vector extensions via `src/utils/simd.h` — auto-scales between SSE (128-bit), AVX (256-bit), and AVX-512 (512-bit) through `VEC_BYTES`.
- **Profiling:** Compile with `-DIHSNPEAKS_PROFILE=1` to enable cycle-accurate phase profiling. Thread-local accumulators get flushed atomically. The `profile_report()` function dumps timing and working-set estimates to stderr.
- **Memory management:** Uses system allocator by default. Optionally links mimalloc (`MIMALLOC=1`). For simplicitly within the hot path, the code allocates and reuses pre-sized buffers per thread rather than frequent malloc/free.
- **Batch mode:** When a directory is provided as target, files are inspected to size the NuFFT plan based on the largest file. Parallel batch processing uses a kthread thread pool.
- **Testing:** No test framework is currently set up in the repository. The `src/debug/` directory contains Python debugging/viewer scripts.
- **No warnings policy:** Compilation with `-Wall -Wextra` is expected to be clean. The makefile adds specific `-Wno-*` flags for certain compiler quirks.
- **Version:** Defined as `IHSNPEAKS_VERSION "v1.1.0-preview"` (work in progress) in `src/main.c`.
