# ihsnpeaks

`ihsnpeaks` is a C command-line periodogram utility for astronomical light-curve analysis. The 1.0 release line is focused on the CLI: high-performance multi-harmonic IHS/Rayleigh-style periodograms, AoVMH(W)/FastChi^2 periodogram's implementation, and Gaussian blur-based Supersmoother-like periodogram, with ANCOVA and F-test for inequality of variance available for reevaluation, as well, as batch processing of OGLE-format `.dat` photometric files (with native multithreading support).

The tool is meant first and foremost as a more performant (and stable) replacement for [FNPEAKS/MPEAKS](http://helas.astro.uni.wroc.pl/deliverables.php?lang=en&active=fnpeaks), [aovdist](https://users.camk.edu.pl/alex/#software), and [FastChi2 v1.03](https://web.archive.org/web/20250525051459/http://public.lanl.gov/palmer/fastchi.html) CLI utilities.

The repository additionally contains `photview` and `spec_viewer` Python scripts meant for debugging purposes, which are not a part of the release

## Usage

```sh
./ihsnpeaks target fmax [options]
```

Self-contained example usec ase:
```sh
wget https://www.astrouw.edu.pl/ogle/ogle4/OCVS/BLAP/phot/phot_ogle4/I/OGLE-BLAP-035.dat
./ihsnpeaks test_data/OGLE-BLAP-035.dat 1000 -n3
```


## Features

- Built-in (NU)FFT backend; FFTW3 is no longer downloaded or linked.
- Full static musl Linux x86-64 and ARM64 release builds with runtime dispatch, plus a Docker-built macOS ARM64 artifact.
- Runtime dispatch variants for `x86-64`, `x86-64-v2`, `x86-64-v2+avx`, `x86-64-v3`, and `x86-64-v4` when supported by the build host.
- ARM64 release dispatch variants for generic ARM64, NEON, SVE 128/256/512, and SVE2 128/256/512.
- Native source builds prefer GNU C23, with fallbacks to GNU C11 and GNU C99 standards.
- C99-compatible aligned allocation fallback is kept for systems without C11 `aligned_alloc`.

## Building from source

The default source build targets the current machine:

```sh
git clone --depth=1 https://github.com/kkrutkowski/ihsnpeaks
cd ihsnpeaks
make native
```
Which can be 'installed' into the system path by following it with
```sh
sudo make install
```

`make` currently is an alias for `make native`. The native build requires a POSIX-like environment, `make`, `readlink`, and a GNU-compatible C compiler. Newer compilers are expected to produce more performant code, but fallbacks from the default `GNU23` C standard to `GNU11` and `GNU99` are included in the codebase..


The release-type multiple dispatch binaries can be compiled by:
```sh
make release-full
```
Which uses Docker to save multiple dispatch `ihsnpeaks-linux-x86_64`, `ihsnpeaks-linux-arm64` and `ihsnpeaks-macos` into `dist/` subdirectory.

This writes the two Linux artifacts plus the universal macOS binary `dist/ihsnpeaks-macos`. Use `make release-macos` to build only the macOS artifact. The macOS build uses `zig cc -target x86_64-macos` and `zig cc -target aarch64-macos`, links project code directly into each thin executable, and merges them with `llvm-lipo`. The only normal runtime dependency is macOS `libSystem`. `MACOS_MIN_VERSION` defaults to `12.0`; set `MACOS_SDK_PATH` when an explicit SDK sysroot is required.

Useful options:

- `-d, --degree, --terms`: number of harmonics.
- `-m, --mode`: peak evaluation/refinement mode.
- `-g, --grid`: periodogram method: default `ihs` or phase-coherent `aov`, `aovmh`, `aobmhw`, `chi`, `chi2`, and `fastchi2` (all phase-coherent keys are aliases of the same method).
- `-s, --save, --spectrum`: write generated spectra to `.tsv` files.
- `-j, --jobs`: limit worker threads (default:0 - unlimited threads)

Run `./ihsnpeaks --help` for the full option list.

## Release Notes

- The 1.0 CLI no longer depends on FFTW3 or mimalloc by default.
- The Linux release binaries use musl static linking and runtime dispatch for x86-64 and ARM64.
- The macOS release binary is a universal Mach-O executable linked through the macOS `libSystem` ABI.
- Dispatch can be forced for testing with `IHSNPEAKS_DISPATCH=<variant>`.
- Stale local binaries are not release artifacts; rebuild with `make native` or `make release`.

## Credits

1. [astropy](https://github.com/astropy/astropy)
2. [aovdist](https://users.camk.edu.pl/alex/)
3. [FastChi2](https://web.archive.org/web/20260000000000*/https://public.lanl.gov/palmer/fastchi.html)
4. [fnpeaks](http://helas.astro.uni.wroc.pl/deliverables.php?active=fnpeaks)
5. [klib](https://github.com/attractivechaos/klib)
6. [FastTransforms.jl](https://github.com/JuliaApproximation/FastTransforms.jl)
7. [qfits](https://www.eso.org/sci/software/eclipse/qfits/html/index.html)
8. [GNU Scientific Library](https://www.gnu.org/software/gsl/)
9. [fast_convert](https://github.com/hermantb/fast_convert)
