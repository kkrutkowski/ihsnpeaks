# ihsnpeaks

`ihsnpeaks` is a C command-line periodogram utility for astronomical light-curve analysis. The 1.0 release line is focused on the CLI: high-performance multi-harmonic IHS/Rayleigh-style periodograms, AoVMH(W)/FastChi^2 periodogram's implementation, and Gaussian blur-based Supersmoother-like periodogram, with ANCOVA and F-test for inequality of variance available for reevaluation, as well, as batch processing of OGLE-format `.dat` photometric files (with native multithreading support).

The tool is meant first and foremost as a more performant (and stable) replacement for [FNPEAKS/MPEAKS](http://helas.astro.uni.wroc.pl/deliverables.php?lang=en&active=fnpeaks), [aovdist](https://users.camk.edu.pl/alex/#software) and [FastChi2 v1.03](https://web.archive.org/web/20250525051459/http://public.lanl.gov/palmer/fastchi.html) CLI utilities.

The repository additionally contain `photview` and `spec_viewer` Python utilities meant for debugging purposes, which are not a part of 

## Features

- Built-in PSWF NuFFT backend; FFTW3 is no longer downloaded or linked.
- No mimalloc dependency in the default build path.
- Full static musl Linux x86-64 release build with runtime dispatch.
- Runtime dispatch variants for `x86-64`, `x86-64-v2`, `x86-64-v2+avx`, `x86-64-v3`, and `x86-64-v4` when supported by the build host.
- Native source builds prefer GNU C23, with fallback to GNU C11 and GNU C99 standards.
- C99-compatible aligned allocation fallback is kept for systems without C11 `aligned_alloc`.

## Building

The default source build targets the current machine:

```sh
make native
```

`make` is an alias for `make native`. The native build requires a POSIX-like environment, `make`, `readlink`, and a GNU-compatible C compiler. Newer compilers generally produce faster binaries, but the makefile probes the compiler and chooses `-std=gnu23`, `-std=gnu11`, or `-std=gnu99`.

Install the native binary with:

```sh
sudo make install
```

The static Linux x86-64 release binary is generated with Docker:

```sh
make release
```

This writes `dist/ihsnpeaks-linux-x86_64`. `dist/` is ignored by git, so release binaries should be regenerated as part of the release process rather than treated as tracked source artifacts.

## Usage

```sh
./ihsnpeaks target fmax [options]
```

Examples:

```sh
./ihsnpeaks test_data/OGLE-BLAP-035.dat 10 -n 3 -d 5
./ihsnpeaks /path/to/photometry_directory 40 -d 1 -m 0 -j 16
./ihsnpeaks target.dat 20 --save --degree 10
```

Useful options:

- `-d, --degree, --terms`: number of harmonics.
- `-m, --mode`: peak evaluation/refinement mode.
- `-g, --grid`: periodogram method: default `ihs` or phase-coherent `aov`, `aovmh`, `aobmhw`, `chi`, `chi2`, and `fastchi2` (all phase-coherent keys are aliases of the same method).
- `-s, --save, --spectrum`: write generated spectra to `.tsv` files.
- `-j, --jobs`: limit worker threads (default:0 - unlimited threads)

Run `./ihsnpeaks --help` for the full option list.

## Release Notes

- The 1.0 CLI no longer depends on FFTW3 or mimalloc by default.
- The release binary uses musl static linking and x86-64 runtime dispatch.
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
