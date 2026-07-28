# Radix-4 / Radix-8 Feasibility Report for nufft1.c

Deep-research analysis (101 agents, 19 sources, 25 claims verified) of whether adding
radix-4 (and radix-8) butterfly steps to the Sande-Tukey DIF FFT in `src/nufft/nufft1.c`
is worthwhile.

---

## 1. Scope of changes for radix-4

**Significantly more than just splitting the butterfly loop.** Three subsystems need rework:

- **Twiddle layout** — `generate_fft_buffer` produces one twiddle per butterfly (W^k).
  Radix-4 needs three (W^k, W^2k, W^3k) with step decreasing 4x per stage.
- **Intra-vector pass** — the precomputed `permutations`/`inv_permutations` shuffle masks
  and `real_twiddles`/`imag_twiddles` arrays are sized for log2(VECF_LEN) radix-2 stages.
  A radix-4 intra-vector pass needs 4-way shuffle masks and a DFT-4 core.
- **Mixed-radix branching** — pure radix-4 requires log2(N) even; odd log2(N) needs one
  leftover radix-2 stage, adding control flow.

The arithmetic savings are real but modest: 25% fewer complex multiplies, ~15% total
operation reduction. The **dominant practical benefit is halving memory passes**, which
matters mainly when the FFT exceeds cache size.

## 2. Bit-reverse permutation

**Cannot remain unchanged.** Radix-4 DIF naturally produces **base-4 digit-reversed**
output, which is provably different from binary bit-reversal. For N=16, they differ on
12 of 16 indices. Mathematically: bit-reversal reverses individual bits
(b3,b2,b1,b0)->(b0,b1,b2,b3), while digit-reversal reverses 2-bit groups as units
(b3b2,b1b0)->(b1b0,b3b2). Both `table_shuffle` and `cobra_shuffle` would need redesign.

One mitigation exists: modifying the radix-4 butterfly (swapping B/C outputs) to emit
bit-reversed order directly, at the cost of extra shuffles inside the butterfly.

## 3. Radix-8: worth it?

**Probably not for this codebase.** A GPU benchmark showed ~22% additional speedup of
radix-8 over radix-4, but that was on Apple Silicon with 128 GPRs/thread. On x86 with
AVX2 (16 YMM registers), a radix-8 butterfly needs 8 vector pairs (16 registers) for
data alone, leaving zero scratch — register spilling would likely negate gains. FFTW
achieves radix-32 codelets via machine code generation, not hand-written SIMD.

More critically: **in the NUFFT pipeline the FFT is only ~1/3 of total runtime**
(interpolation/spreading dominates at ~2x the FFT cost). Even an optimistic 40% FFT
speedup yields at most ~8-13% end-to-end improvement — a modest return given the
implementation complexity.

## Alternative approaches worth considering

- **Stockham autosort (out-of-place) radix-2** — eliminates the explicit permutation
  pass entirely without changing radices. Potentially a simpler optimization path than
  switching to radix-4.
- **Split-radix** (one N/2 + two N/4 sub-DFTs) — achieves the theoretical minimum
  arithmetic count while preserving standard bit-reversal ordering and requiring less
  restructuring of the existing code.

## Caveats

1. Register pressure data (radix-8 at 38 registers) comes from Apple Silicon GPU
   (128 GPRs/thread), not x86 CPU SIMD. On x86 with AVX2/AVX-512 the landscape is
   fundamentally different and likely more constrained.
2. The 22% radix-8-over-radix-4 speedup benchmark is from a GPU kernel (arxiv preprint),
   not a CPU SIMD implementation using GCC vector extensions.
3. The NUFFT timing ratio (interpolation ~2x FFT) is from a 2003 paper; modern hardware
   may shift this ratio, though recent GPU measurements (2025) suggest interpolation
   dominance has if anything increased.
4. No benchmark data exists for radix-4 specifically within this codebase's architecture
   (GCC vector extensions, interleaved real/imag arrays, intra-vector passes).
5. The FFTW observation that minimal-arithmetic implementations are suboptimal applies to
   general-purpose libraries, not necessarily to this specialized single-size-per-plan
   NUFFT context.

## Open questions

- What FFT sizes are actually used in practice? If consistently small (<=1024) and
  fitting in L1, the memory-pass-reduction benefit of radix-4 is largely negated.
- Would Stockham autosort provide comparable or better speedup than radix-4, given that
  the current cobra_shuffle/table_shuffle is a separate cache-unfriendly pass?
- For the specific VEC_BYTES configuration, how many SIMD registers does the current
  radix-2 butterfly consume, and would radix-4 cause spilling?
- Could split-radix be a better fit than pure radix-4?

## Sources

- [TI SPRA152 — Radix-4 FFT](https://www.ti.com/lit/an/spra152/spra152.pdf)
- [NXP AN3666 — FFT implementation](https://www.nxp.com.cn/docs/en/application-note/AN3666.pdf)
- [Fabian Giesen — Notes on FFTs for implementers (2023)](https://fgiesen.wordpress.com/2023/03/19/notes-on-ffts-for-implementer/)
- [FFTW paper (IEEE)](https://www.fftw.org/fftw-paper-ieee.pdf)
- [FFTW split-radix paper](https://www.fftw.org/newsplit.pdf)
- [US Patent 5473556A — Digit reversal](https://patents.google.com/patent/US5473556A/en)
- [Fessler & Sutton 2003 — NUFFT](https://web.eecs.umich.edu/~fessler/papers/files/jour/03/pre/tsp,nufft.pdf)
- [arxiv 2603.27569 — GPU radix-8 benchmark](https://arxiv.org/html/2603.27569v1)
- [arxiv 2605.10678 — GPU NUFFT timing (2025)](https://arxiv.org/html/2605.10678)
- [DSP StackExchange — Bit reversal necessity](https://dsp.stackexchange.com/questions/5938/why-exactly-is-a-bit-radix-reverse-required-when-calculating-the-fft)
- [DSP StackExchange — Higher radix FFTs](https://dsp.stackexchange.com/questions/79114/why-not-use-higher-power-of-2-radix-fft)

---

## Appendix: COBRAVO Implementation Note

Reference: Lokhmotov & Mycroft, "Brief Announcement: Optimal Bit-Reversal Using Vector
Permutations," SPAA '07. See `https://dl.acm.org/doi/pdf/10.1145/1248377.1248411`.

### What COBRAVO is

COBRAVO = COBRA (cache-optimal tiling) + BRAVO (intra-tile bit-reversal via SIMD vector
interleave permutations). It replaces the **scalar element-by-element swaps** inside each
COBRA tile with hardware vector permute instructions. The bit-reverse lookup table is
retained for tile-level index computation and for the final vector reordering step.

### Mapping to existing code

The current `cobra_shuffle` in `nufft1.c` already implements the COBRA tiling:

- `LOG_BLOCK_WIDTH = 6`, `BLOCK_WIDTH = 64` (the tile dimension B)
- Buffer: `cobra_real` / `cobra_imag`, each `BLOCK_WIDTH * BLOCK_WIDTH = 4096` floats
- Outer loop iterates over the middle field `b` (the `b_size` iterations)
- For each `b`: loads a 64x64 tile into the buffer, performs scalar swaps, writes back

**What changes:** only the middle section — the scalar swap loops (the three nested loops
that do `float tr = real[v_idx]; real[v_idx] = buffer_real[b_idx]; ...`) are replaced
with the BRAVO vector permutation sequence described below. Everything else (tile loading,
tile writing, tile-level index computation via `reverse_bits`) stays the same.

### The BRAVO intra-tile permutation algorithm

Given a tile of B^2 elements stored contiguously in `T[B*B]` (row-major), viewed as
`N_vec = B*B / W` SIMD vectors of width `W = VECF_LEN`:

**Step 1: Interleave stages**

Perform `log2(N_vec)` stages. At stage `s` (s = 0, 1, ..., log2(N_vec)-1):

```
stride = 1 << s
for group_start = 0; group_start < N_vec; group_start += 2 * stride:
    for i = group_start; i < group_start + stride; i++:
        A = vec[i]
        B = vec[i + stride]
        vec[i]          = interleave_low(A, B)
        vec[i + stride] = interleave_high(A, B)
```

Where:
- `interleave_low(A, B)`  = {A[0], B[0], A[1], B[1], ..., A[W/2-1], B[W/2-1]}
- `interleave_high(A, B)` = {A[W/2], B[W/2], A[W/2+1], B[W/2+1], ..., A[W-1], B[W-1]}

**Step 2: Bit-reverse the vector ordering**

After all interleave stages, the *contents* of each vector are correct, but the vectors
themselves are in the wrong order. Permute the N_vec vectors by bit-reversing their index:

```
for i = 0; i < N_vec; i++:
    j = reverse_bits(i, log2(N_vec))
    if j > i:
        swap vec[i] and vec[j]
```

This uses the existing `bit_reverse_table` (or the `reverse_bits` helper).

### Concrete example (N=16, W=4)

Input vectors: V0=[0,1,2,3], V1=[4,5,6,7], V2=[8,9,10,11], V3=[12,13,14,15]
N_vec = 4, log2(N_vec) = 2 stages.

Stage s=0 (stride=1): pair (0,1) and (2,3):
  V0, V1 = il(V0,V1), ih(V0,V1) = [0,4,1,5], [2,6,3,7]
  V2, V3 = il(V2,V3), ih(V2,V3) = [8,12,9,13], [10,14,11,15]

Stage s=1 (stride=2): pair (0,2) and (1,3):
  V0, V2 = il(V0,V2), ih(V0,V2) = [0,8,4,12], [1,9,5,13]
  V1, V3 = il(V1,V3), ih(V1,V3) = [2,10,6,14], [3,11,7,15]

After interleaving: [0,8,4,12], [1,9,5,13], [2,10,6,14], [3,11,7,15]

Bit-reverse vector order (2-bit reversal: 0->0, 1->2, 2->1, 3->3):
  Position 0 <- vector 0: [0,8,4,12]
  Position 1 <- vector 2: [2,10,6,14]
  Position 2 <- vector 1: [1,9,5,13]
  Position 3 <- vector 3: [3,11,7,15]

Final output: [0,8,4,12, 2,10,6,14, 1,9,5,13, 3,11,7,15]
Verification: rev4(0)=0, rev4(1)=8, rev4(2)=4, rev4(3)=12, ... correct.

### Expressing interleaves with NANOFFT_SHUFFLE2

The existing `NANOFFT_SHUFFLE2(a, b, mask)` macro maps directly to interleave operations.
The mask indices 0..W-1 select from `a`, indices W..2W-1 select from `b`.

For W = VECF_LEN = 8 (AVX2, 256-bit):
```c
// interleave_low: take elements 0..3 from each, interleaved
#define IL_LOW_MASK  {0, 8, 1, 9, 2, 10, 3, 11}
// interleave_high: take elements 4..7 from each, interleaved
#define IL_HIGH_MASK {4, 12, 5, 13, 6, 14, 7, 15}

VECF lo = NANOFFT_SHUFFLE2(A, B, (VECF_INT)IL_LOW_MASK);
VECF hi = NANOFFT_SHUFFLE2(A, B, (VECF_INT)IL_HIGH_MASK);
```

For W = 16 (AVX-512, 512-bit):
```c
#define IL_LOW_MASK  {0,16, 1,17, 2,18, 3,19, 4,20, 5,21, 6,22, 7,23}
#define IL_HIGH_MASK {8,24, 9,25, 10,26, 11,27, 12,28, 13,29, 14,30, 15,31}
```

For W = 4 (SSE, 128-bit):
```c
#define IL_LOW_MASK  {0, 4, 1, 5}
#define IL_HIGH_MASK {2, 6, 3, 7}
```

### Parameters for the existing code

With `BLOCK_WIDTH = 64` and `VECF_LEN = W`:

| Parameter       | Value (W=8) | Value (W=16) |
|-----------------|-------------|--------------|
| B (tile dim)    | 64          | 64           |
| B^2 (tile size) | 4096        | 4096         |
| N_vec           | 512         | 256          |
| Interleave stages | 9         | 8            |
| Pairs per stage | 256         | 128          |
| Vector swaps (step 2) | ~256 | ~128         |

### Operation counts vs current scalar approach

Current scalar `cobra_shuffle` inner loops: up to B^2/2 = 2048 conditional swaps, each
swap = 2 loads + 2 stores (for real) + 2 loads + 2 stores (for imag) = 8 scalar memory
operations. Total: ~16384 scalar memory ops per tile.

COBRAVO per tile (W=8):
- 9 stages x 256 pairs x 2 shuffles = 4608 vector shuffle instructions (each processes
  8 floats)
- ~256 vector swaps (each = 2 vector loads + 2 vector stores)
- Total: ~5120 vector instructions processing 8 floats each

The speedup comes from: (1) each shuffle instruction replaces multiple scalar
load/store/compare operations, (2) no data-dependent branching (all pairs are processed
unconditionally), (3) better instruction-level parallelism since shuffles are pure
register operations.

### Pseudocode for the replacement function

```c
static inline void bravo_shuffle_tile(float *tile_real, float *tile_imag, int n_vec) {
    // tile_real/tile_imag point to B*B contiguous floats (the buffer T)
    // n_vec = B*B / VECF_LEN
    int log_nvec = intlog2((uint32_t)n_vec);

    // Step 1: interleave stages
    for (int s = 0; s < log_nvec; ++s) {
        int stride = 1 << s;
        for (int group = 0; group < n_vec; group += 2 * stride) {
            for (int i = group; i < group + stride; ++i) {
                int j = i + stride;
                VECF ar = LOAD_VEC(&tile_real[i * VECF_LEN]);
                VECF br = LOAD_VEC(&tile_real[j * VECF_LEN]);
                VECF ai = LOAD_VEC(&tile_imag[i * VECF_LEN]);
                VECF bi = LOAD_VEC(&tile_imag[j * VECF_LEN]);

                STORE_VEC(&tile_real[i * VECF_LEN], NANOFFT_SHUFFLE2(ar, br, IL_LOW_MASK));
                STORE_VEC(&tile_real[j * VECF_LEN], NANOFFT_SHUFFLE2(ar, br, IL_HIGH_MASK));
                STORE_VEC(&tile_imag[i * VECF_LEN], NANOFFT_SHUFFLE2(ai, bi, IL_LOW_MASK));
                STORE_VEC(&tile_imag[j * VECF_LEN], NANOFFT_SHUFFLE2(ai, bi, IL_HIGH_MASK));
            }
        }
    }

    // Step 2: bit-reverse the vector ordering
    for (int i = 0; i < n_vec; ++i) {
        int j = (int)reverse_bits((uint32_t)i, (uint32_t)log_nvec);
        if (j > i) {
            // swap vectors i and j (both real and imag)
            VECF tmp_r = LOAD_VEC(&tile_real[i * VECF_LEN]);
            STORE_VEC(&tile_real[i * VECF_LEN], LOAD_VEC(&tile_real[j * VECF_LEN]));
            STORE_VEC(&tile_real[j * VECF_LEN], tmp_r);
            VECF tmp_i = LOAD_VEC(&tile_imag[i * VECF_LEN]);
            STORE_VEC(&tile_imag[i * VECF_LEN], LOAD_VEC(&tile_imag[j * VECF_LEN]));
            STORE_VEC(&tile_imag[j * VECF_LEN], tmp_i);
        }
    }
}
```

### Integration into cobra_shuffle

Replace the three scalar-swap loop sections in `cobra_shuffle` (the sections guarded by
`if (less) { ... swap ... }`) with a single call to `bravo_shuffle_tile` operating on the
buffer after it has been loaded. The tile loading (first loop, copying from `real[]`/
`imag[]` into `buffer_real[]`/`buffer_imag[]`) and tile writing (third loop, copying back
to the main arrays at the transposed tile location) remain unchanged.

The `bit_reverse_table` continues to be used for:
1. Computing `b_rev = reverse_bits(b, num_b_bits)` (tile-level index mapping)
2. The vector reordering in Step 2 of BRAVO (via `reverse_bits`)

### Note on radix-4 interaction

If radix-4 stages are added to the FFT butterfly (producing base-4 digit-reversed output),
the COBRAVO tiling structure remains valid — it is agnostic to which permutation is being
computed. However, the BRAVO interleave sequence implements binary bit-reversal
specifically. For digit-reversal, either:

(a) Use the "modified butterfly" approach (swap B/C outputs in the radix-4 core) so the
    FFT still emits bit-reversed order, and COBRAVO works unchanged, OR
(b) Replace the BRAVO interleave sequence with one that computes digit-reversal. The
    structure is identical (same stages, same strides), but the shuffle masks differ:
    instead of interleaving individual elements, you interleave 2-element groups (pairs),
    preserving intra-pair order. For W=8 this means:
    `IL_LOW_MASK_DIGIT = {0, 1, 8, 9, 2, 3, 10, 11}` (interleave pairs, not singles).
