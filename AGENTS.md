# ihsnpeaks — Project Guide

## Project Overview

`ihsnpeaks` is a high-performance C command-line periodogram utility for astronomical light-curve analysis, licensed under GPLv3. It is designed as a drop-in replacement for legacy tools (FNPEAKS/MPEAKS, aovdist, FastChi2).

**Key capabilities:**
- IHS/Rayleigh-style multi-harmonic periodograms
- AoV(MH), AoVMH(W), FastChi² periodogram implementations
- Gaussian blur-based Supersmoother-like periodogram (gbls/gbaw)
- BLS (Box Least Squares) transit search
- Bayesian statistic (Relative Evidence Ratio) and Classical Raw NLL evaluation with Pochhammer rising factorial series expansions
- ANCOVA and F-test for inequality-of-variance reevaluation
- Batch processing of OGLE-format `.dat` photometric files
- MAST/TESS FITS light-curve input (PDCSAP_FLUX/SAP_FLUX/FLUX via qfits), mixed `.dat`+`.fits` batches
- hwloc-based L3-cache-aware thread pools with per-PU CPU binding and first-touch NUMA placement
- Native multithreading (pthreads) with L3-aware per-pool chunked dispatch and work-stealing
- Built-in PSWF-based NuFFT backend (no FFTW3 dependency)

**Architecture:**
- Single-file CLI binary: `src/main.c` is the entry point
- Core periodogram logic lives in `src/process.h` (inline header)
- NuFFT backend in `src/nufft/nufft1.c` and `nufft1.h`
- SIMD vectorization via GCC/Clang vector extensions (`src/utils/simd.h`)
- Uses klib (kvec, kthread, ketopt, kopen) for data structures and helpers
- Uses sds (Redis dynamic strings) for string handling
- Uses qfits (header-only) for FITS light-curve reading (`src/utils/readout.h` `read_fits()`)
- Uses hwloc (statically linked, built from source) for L3 topology discovery and CPU binding (`src/hwloc_topo.h`)
- L3-aware dispatcher (`src/pool.h`) creates one `kt_forpool` per L3 die, distributes files in chunks via a shared atomic counter
- Runtime dispatch: single binary selects the optimal CPU feature variant at startup
- Static musl builds for Linux; universal Mach-O binary for macOS

**Repository structure:**

```
.
├── src/
│   ├── main.c           — Entry point, CLI argument parsing, batch orchestration
│   ├── params.h         — Parameter struct, init, parsing, help text
│   ├── metadata.h       — File discovery, NuFFT plan sizing, plan cache init
│   ├── process.h        — Core periodogram evaluation pipeline (main logic)
│   ├── profile.h        — Instrumentation/profiling macros (IHSNPEAKS_PROFILE)
│   ├── hwloc_topo.h     — hwloc topology probing, L3 domain enumeration, per-PU CPU binding, XML export to ~/.ihsnpeaks/hwloc_config
│   ├── pool.h           — L3-aware thread-pool dispatcher: per-pool buffer sets, chunked atomic file dispatch
│   ├── nufft/           — PSWF-based non-uniform FFT backend
│   │   ├── nufft1.c     — NuFFT implementation
│   │   ├── nufft1.h     — NuFFT public API
│   │   ├── scaling.c    — PSWF scaling coefficients generator
│   │   └── scaling_generic.h — Fallback PSWF scaling coefficients
│   ├── utils/           — Shared utilities
│   │   ├── common.h     — Core types: buffer_t, parameters, peak_t, enums
│   │   ├── simd.h       — GCC vector extension SIMD wrappers (ln_ps, correctPower, rising factorials)
│   │   ├── aov.h        — AoV-based periodogram helpers
│   │   ├── convolution.h — Gaussian blur / BLS direct evaluation
│   │   ├── readout.h    — Buffer allocation, I/O (read_dat/read_fits), peak output formatting
│   │   ├── compat.h     — C11 aligned_alloc fallback for C99
│   │   ├── trig.h       — Trig utilities (sin2pif_tls, cos2pif_tls)
│   │   ├── spectrum.h   — Spectrum column capture
│   │   └── strings.h    — String formatting helpers
│   └── debug/           — Debugging helpers (not part of release)
│       ├── photview     — Python photometry viewer
│       └── spec_viewer.py — Python spectrum viewer
├── include/
│   ├── klib/            — External klib headers (kvec, kthread, ketopt, kopen)
│   ├── qfits/           — qfits library header (header-only, reads MAST FITS tables)
│   ├── fast_convert.h   — Fast float parsing
│   ├── fdist.h          — F-distribution CDF & Pochhammer series
│   ├── json.h           — sheredom JSON parser
│   └── sds.h            — Redis SDS dynamic strings
├── downstream/
│   └── lc-qt/           — Qt6 GUI light curve analysis application
├── dispatch/
│   ├── variants.py      — Runtime dispatch variant definitions
│   ├── build_release.py — x86-64 musl release builder (builds hwloc per target)
│   ├── build_release_arm64.py — ARM64 musl release builder (builds hwloc per target)
│   ├── build_release_macos.py — macOS universal binary builder (per-arch hwloc)
│   └── zig_musl_cc.sh   — Zig-based cross-compilation wrapper
├── Dockerfile.release   — Multi-stage Docker build for release artifacts
├── makefile             — Primary build system (GNU Make; fetches+builds hwloc from source)
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

**hwloc is built from source (v2.14.0) for each release target.** Each `dispatch/build_release*.py` downloads the hwloc tarball into the variant's build directory and compiles `libhwloc.a` with the target cross-compiler (zig musl for Linux, zig darwin for macOS per-arch slices). The native makefile follows the same pattern in `build/native/<allocator>/`.

### Profile build

```sh
make profile         # Builds a profiled binary with -DIHSNPEAKS_PROFILE=1
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

- `target`: path to a `.dat` or `.fits` file, or a directory of `.dat`/`.fits` files (mixed batches are supported)
- `fmax`: upper frequency bound for the periodogram grid

**Common options:**

| Flag | Description |
|------|-------------|
| `-d N` | Number of harmonics/terms (default: 3) |
| `-m N` | Peak evaluation mode 0–6 (default: 2) |
| `-g METHOD` | Periodogram method: `ihs` (default for `-d <= 4`), `aov`/`aovmh`/`chi`/`chi2` (default for `-d >= 5`) |
| `--statistic TYPE` | Statistic type: `bayes` (Relative Evidence Ratio) or `raw` (Classical NLL; default: `bayes` for m0-m2 AoV, `raw` otherwise) |
| `-e METHOD` | Peak evaluation: `gbls[alpha]`, `gbaw[alpha]`, or `bls[a,b,count]` |
| `-t THRESHOLD` | Detection threshold in dB (default: 10.0) |
| `-n N` | Max number of peaks to report (default: 10) |
| `-s` | Save full spectrum to `.tsv` files |
| `-p` | Prewhiten (subtract detected signals before finding next peaks) |
| `-j N` | Limit worker threads (default: 0 = all available PUs) |
| `-b MODE` | Thread binding: `0`/`false`, `1`/`strict`, `2`/`auto` (default), `3`/`cache` |
| `-o N` | Oversampling factor (default: 5.0) |
| `-f FMIN` | Lower frequency bound (default: 0.0) |
| `--period` | Output periods instead of frequencies |
| `--nufft MODE` | NuFFT backend: `43`/`pswf43` (default) or `21`/`pswf21` |
| `--generate` | Probe hardware topology via hwloc, save to `~/.ihsnpeaks/hwloc_config`, and exit |
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
- **Memory management:** Uses system allocator by default. Optionally links mimalloc (`MIMALLOC=1`). For simplicity within the hot path, the code allocates and reuses pre-sized buffers per thread rather than frequent malloc/free.
- **Batch mode:** When a directory is provided as target, the metadata pass discovers all matching `.dat`/`.fits` files, but **inspects only the largest `.dat` file and the largest `.fits` file** (one per format, by `stat` size) to size the NuFFT plan. Buffer sizes are derived from these two files only.
- **Batch metadata scan performance (IMPORTANT):** The metadata pass must NOT loop over all files calling `inspect_dat_file()` or `inspect_fits_file()`. Doing so introduces severe startup latency on large directories (e.g., 10+ seconds on 15 000+ files where < 1 second is acceptable). The current implementation in `metadata.h::process_path()` tracks `largest_dat_path` / `largest_fits_path` during `readdir()` and opens at most two files after the loop. Any refactor must preserve this property.
- **Threading / hwloc:** At startup, `load_or_probe_topology()` probes the hardware via hwloc (or loads the cached topology from `~/.ihsnpeaks/hwloc_config`). One worker thread pool is created per active L3 domain; workers are pinned to specific PUs via `hwloc_set_cpubind(HWLOC_CPUBIND_THREAD)` before their per-thread `buffer_t` is allocated, so first-touch NUMA placement puts the per-thread buffers on the local node. Shared NuFFT twiddle arrays (`nufftTwiddleReal/Imag`) are allocated before workers exist, so the OS kernel distributes them across NUMA nodes by default (acceptable for read-mostly data).
- **Testing:** No test framework is currently set up in the repository. The `src/debug/` directory contains Python debugging/viewer scripts.
- **No warnings policy:** Compilation with `-Wall -Wextra` is expected to be clean. The makefile adds specific `-Wno-*` flags for certain compiler quirks.
- **Version:** Defined as `IHSNPEAKS_VERSION "v1.1.0-preview"` (work in progress) in `src/main.c`.
