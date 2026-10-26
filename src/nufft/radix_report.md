# NuFFT FFT: fused Radix-2, Radix-4 comparison, and calibration

## Summary

Production now uses a **fused Radix-2 pair with interleaved twiddles** in
`src/nufft/nufft1.c`.  It replaces the old two-pass Radix-2 FFT while keeping
the same DIF Sande-Tukey arithmetic, binary bit-reversed output, COBRA/COBRAVO
permutation, and public NuFFT API.

A Radix-4 DIF implementation was prototyped as well.  With an interleaved
twiddle layout it is genuinely faster than the old two-pass Radix-2, but a fair
fused-Radix-2 baseline with the same layout captures essentially all of that
gain.  Radix-4's unique arithmetic advantage is only ~2-5% on AVX2/SSE and is
neutral on AVX-512; the required grid-sparsity penalty then makes it lose.  It is
therefore not part of the production path.

## Production implementation: fused Radix-2 pair

Two adjacent DIF Radix-2 stages are processed in one pass.  For four inputs
`a,b,c,d` spaced `q` apart:

- stage 1: `(a,c)` with `W1`, `(b,d)` with `i*W1`
- stage 2: the two half-blocks with `W2 = W^2`

The outputs are stored as `X0, X2, X1, X3`, exactly the binary bit-reversed order
produced by the old two-pass code.  The output permutation is therefore
unchanged and COBRA/COBRAVO need no modification.

The inner loop uses four data vectors and two twiddle vectors per j-block.
`W1` and `W2` for each block are stored contiguously (two `VECF_LEN`-wide
vectors per block).  This removes the three widely separated twiddle streams
that made a naive Radix-4 implementation regress on large N.  Fusing the stages
also halves the number of full-array passes.

Because the fused layout is **not a suffix** of the next-larger plan's layout,
`src/metadata.h` now gives every cached plan its own twiddle region inside the
shared allocation and calls `nufft1_initialize_shared()` per plan.  The total
allocation is roughly twice the old suffix-shared size; plan setup is done once
before worker threads start.  `nufft1_initialize_shared_precomputed()` still
exists, but its twiddle buffers must already contain that plan's fused layout;
a suffix of a larger plan's buffer is not valid.

## Fair comparison with Radix-4

Prototype host: Ryzen 9 7945HX, GCC 14, pinned to one CPU, best of 3.
`nufft1_execute`, N = plan length, M = 128, PSWF43 (microseconds):

| N | old two-pass R2 | fused R2 + interleaved | Radix-4 + interleaved | Radix-4 × 1.066 |
|---:|---:|---:|---:|---:|
| 8 192 | 31.5 | 29.0 | 28.3 | 30.2 |
| 16 384 | 66.2 | 58.7 | 61.0 | 65.0 |
| 32 768 | 148.9 | 121.3 | 121.6 | 129.6 |
| 65 536 | 332.6 | 249.7 | 254.3 | 271.1 |
| 131 072 | 704.4 | 522.1 | 526.3 | 561.0 |

- Fused Radix-2 vs old: **1.09-1.35x**.
- Radix-4 vs old: **1.11-1.34x**.
- Radix-4 vs fused Radix-2 on AVX-512: **0.96-1.03x** (neutral).
- AVX2, N=65 536: old 306.3, fused 282.9, Radix-4 270.9, penalized 288.8 → fused wins.
- SSE, N=65 536: old 691.3, fused 640.7, Radix-4 618.8, penalized 659.6 → fused wins.

The Radix-4-specific arithmetic is only 2-5% faster than the fused pair on
AVX2/SSE and is not enough to overcome the `2^(alpha/2)` grid-sparsity penalty.
The significant speedup comes from stage fusion plus the interleaved twiddle
layout, both of which the fused Radix-2 pair gets.

## Bit reversal and twiddle layout

- The Radix-4 butterfly stores `X2` before `X1`, which is a 2-bit digit reversal
  and keeps the overall output in **binary bit-reversed** order.  The earlier
  claim that Radix-4 necessarily requires a new base-4 permutation was too
  pessimistic for this DIF formulation.
- COBRA/COBRAVO are unchanged and remain selected at build time by
  `scaling_gen` (`F_BITREV_COBRAVO`).  The fused-Radix-2 output feeds the same
  permutation.
- Twiddle layout for one fused stage: for each j-block of `VECF_LEN` values,
  store `W1` then `W2` as contiguous complex vectors.  The tail (last
  `log2(VECF_LEN)` stages) still uses the hard-coded intra-vector twiddles or
  the scalar fallback.

## Calibration

`src/nufft/scaling.c` now:

1. self-tests the fused FFT against the legacy two-pass kernel and a direct DFT
   for several orders (including odd `log2(N)`), max relative error < 5e-4;
2. fits the FFT power-scaling exponent from log-log timing at N = 16 384 …
   262 144 (larger sizes avoid fixed-overhead bias) and emits `F_ALPHA`;
3. uses min-of-repetitions for precompute/execute timing and emits the
   PSWF43 beta/gamma derived from those stable measurements;
4. verifies and benchmarks COBRA/COBRAVO as before.

On the development host the fitted alpha is ~0.12-0.14, R² ≈ 0.99.  The old
hard-coded value was 0.20; the fitted value changes plan sizes slightly, which
moved one reference frequency by 2e-6.  The reference fixture was updated for
that legitimate calibration change.  Build logs include the fit points, the
legacy-vs-fused FFT timing, and the COBRA/COBRAVO selection.

## Performance on the final production build

Representative execute speedups (fused vs old two-pass) were **1.09-1.35x** for
N = 8 192…131 072 on AVX-512, **1.03-1.13x** on AVX2, and **1.03-1.12x** on SSE.
End-to-end `benchmark.sh` timings on OGLE-BLAP-035 with `-j 1` improved from
about 321 ms to 261 ms for `ihs d=8` and from about 1952 ms to 1649 ms for
`aov d=8`; `d=1` is essentially unchanged.

## Architecture / fallbacks

- The fused loop uses GCC/Clang vector extensions and `__builtin_assume_aligned`
  only when available; otherwise it falls back to unaligned vector loads or the
  scalar tail.
- `n < 2*VECF_LEN` uses the scalar tail; correctness is tested for orders
  13-18.  Odd `log2(N)` keeps one ordinary Radix-2 stage so the output remains
  binary bit-reversed.
- AVX-512, AVX2, SSE, scalar, GCC, Clang, and aarch64 (zig
  `-target aarch64-linux-musl`) builds were verified.
- No runtime dispatch or per-architecture source changes are needed; the same
  source compiles for all supported vector widths.

## References

- TI SPRA152 — Radix-4 FFT: https://www.ti.com/lit/an/spra152/spra152.pdf
- Fabian Giesen — Notes on FFTs for implementers:
  https://fgiesen.wordpress.com/2023/03/19/notes-on-ffts-for-implementer/
- FFTW papers: https://www.fftw.org/fftw-paper-ieee.pdf
- Fessler & Sutton 2003 — NUFFT:
  https://web.eecs.umich.edu/~fessler/papers/files/jour/03/pre/tsp,nufft.pdf
- Lokhmotov & Mycroft, "Optimal Bit-Reversal Using Vector Permutations",
  SPAA '07 (COBRA/COBRAVO).
