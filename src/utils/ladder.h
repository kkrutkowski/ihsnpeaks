#ifndef LADDER_H
#define LADDER_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "../nufft/nufft1.h"
#include "simd.h"
#include "trig.h"

/* =========================================================================
 * Two-level twiddle ladder.
 *
 *  level 0 (state): for every harmonic q = 1..Q, the complex point twiddle of
 *           the current block, already scaled by the per-point amplitude
 *           `mult` (or by `weight = 1/(dy^2 + eps)` when `mult` is NULL).
 *           The state is mutable during a run and is advanced by exactly one
 *           complex multiply per block against the pregenerated per-harmonic
 *           step.
 *  level 1 (anchor): every `anchor` (= MAX_TWIDDLE_REUSE) blocks the state is
 *           regenerated from an exact double-precision phase instead of being
 *           multiplied again.  That bounds the accumulated rounding to the
 *           reuse budget and removes the need for a third level.
 *
 * state_* is written during the run; step_* is generated once, before the
 * block loop, and only read afterwards.  Both are indexed
 * [harmonic * pad + point] with pad a multiple of VEC_LEN.
 *
 * Resync graph: phase = x_i * q * (fmin + df * block_len * block) in double ->
 * DVEC multiply -> fractional part -> one VEC -> sin2pif_vec/cos2pif_vec.
 *
 * MAX_TWIDDLE_REUSE is the level-0 reuse budget; it does not have to be a
 * power of two (the ladder no longer carries a radix decomposition).  Build
 * systems pass -DMAX_TWIDDLE_REUSE=N; the default below only matters for
 * translation units that are built without it.
 * ========================================================================= */

#ifndef MAX_TWIDDLE_REUSE
#    define MAX_TWIDDLE_REUSE 8
#endif

#if MAX_TWIDDLE_REUSE < 2
#    error "MAX_TWIDDLE_REUSE must be >= 2"
#endif

typedef struct {
    uint32_t harmonics;     /* Q: the ladder carries q = 1..Q */
    uint32_t anchor;        /* level-0 reuse budget = MAX_TWIDDLE_REUSE */
    uint32_t block_len;     /* frequencies per block (NuFFT output length) */
    size_t n;               /* number of samples */
    size_t pad;             /* VEC_LEN-aligned n (state/step array stride) */
    double fmin;            /* base frequency of block 0 */
    double df;              /* frequency step */
    const double *x;        /* sample times, n entries */
    const float *mult;      /* per-point amplitude folded into the state; NULL => weight */
    float *weight;          /* pad: 1/(dy^2 + eps), only filled when mult == NULL */
    float *state_r;         /* Q * pad */
    float *state_i;         /* Q * pad */
    float *step_r;          /* Q * pad; NULL when the grid spans a single block */
    float *step_i;          /* Q * pad */
    float *yin_r;           /* pad scratch: y-weighted input of one harmonic */
    float *yin_i;           /* pad scratch */
    uint32_t cap_harmonics; /* allocation bookkeeping */
    size_t cap_pad;
} twiddle_ladder;

static inline void *ladder_aligned_alloc(size_t count, size_t size) {
    if (count == 0U || size == 0U || count > SIZE_MAX / size) return NULL;
    size_t bytes = count * size;
    bytes = (bytes + 63U) & ~(size_t)63U;
    return aligned_alloc(64, bytes);
}

/* Free every ladder array and reset the metadata. */
static inline void ladder_free(twiddle_ladder *ladder) {
    if (!ladder) return;
    free(ladder->state_r);
    free(ladder->state_i);
    free(ladder->step_r);
    free(ladder->step_i);
    free(ladder->weight);
    free(ladder->yin_r);
    free(ladder->yin_i);
    memset(ladder, 0, sizeof(*ladder));
}

/* Grow the arrays to hold `harmonics` rows of `pad` points.  `need_steps`
 * requests the pregenerated per-harmonic step (only needed when the grid spans
 * more than one block).  Arrays already large enough are kept untouched. */
static inline int ladder_alloc(twiddle_ladder *ladder, uint32_t harmonics, size_t pad, bool need_steps) {
    if (!ladder || harmonics == 0U || pad == 0U) return -1;
    size_t count = (size_t)harmonics * pad;
    if (count / (size_t)harmonics != pad) return -1;

    bool rows_ok = ladder->state_r && ladder->state_i && ladder->cap_harmonics >= harmonics && ladder->cap_pad >= pad;
    bool steps_ok = !need_steps || (ladder->step_r && ladder->step_i);
    bool scratch_ok = ladder->weight && ladder->yin_r && ladder->yin_i;

    if (rows_ok && steps_ok && scratch_ok) return 0;

    if (!rows_ok) {
        free(ladder->state_r);
        free(ladder->state_i);
        free(ladder->step_r);
        free(ladder->step_i);
        ladder->state_r = ladder->state_i = NULL;
        ladder->step_r = ladder->step_i = NULL;
        ladder->cap_harmonics = 0;
        ladder->cap_pad = 0;

        ladder->state_r = (float *)ladder_aligned_alloc(count, sizeof(float));
        ladder->state_i = (float *)ladder_aligned_alloc(count, sizeof(float));
        if (!ladder->state_r || !ladder->state_i) {
            ladder_free(ladder);
            return -1;
        }
        ladder->cap_harmonics = harmonics;
        ladder->cap_pad = pad;
    }

    if (need_steps && !(ladder->step_r && ladder->step_i)) {
        ladder->step_r = (float *)ladder_aligned_alloc(count, sizeof(float));
        ladder->step_i = (float *)ladder_aligned_alloc(count, sizeof(float));
        if (!ladder->step_r || !ladder->step_i) {
            ladder_free(ladder);
            return -1;
        }
    }

    if (!ladder->weight) {
        ladder->weight = (float *)ladder_aligned_alloc(pad, sizeof(float));
        if (!ladder->weight) {
            ladder_free(ladder);
            return -1;
        }
    }
    if (!ladder->yin_r || !ladder->yin_i) {
        ladder->yin_r = (float *)ladder_aligned_alloc(pad, sizeof(float));
        ladder->yin_i = (float *)ladder_aligned_alloc(pad, sizeof(float));
        if (!ladder->yin_r || !ladder->yin_i) {
            ladder_free(ladder);
            return -1;
        }
    }
    return 0;
}

/* Exact-phase regeneration of one harmonic state at `block`:
 * state = amp * exp(2*pi*i * x_i * (h+1) * (fmin + df * block_len * block)). */
/* Kept out of line: the vector trig body is large and would otherwise sit in
 * the middle of the per-block sweep loop. */
static void ladder_regenerate_harmonic(twiddle_ladder *ladder, uint32_t h, size_t block) {
    const size_t n = ladder->n;
    const size_t pad = ladder->pad;
    const size_t n_vec = n - (n % (size_t)VEC_LEN);
    const float *restrict amp = ladder->mult ? ladder->mult : ladder->weight;
    const double freq_block = ladder->fmin + ladder->df * ((double)block * (double)ladder->block_len);
    const double scale = (double)(h + 1U) * freq_block;
    float *restrict state_r = ladder->state_r + (size_t)h * pad;
    float *restrict state_i = ladder->state_i + (size_t)h * pad;

    size_t i = 0;
    for (; i < n_vec; i += VEC_LEN) {
        DVEC x_lo, x_hi;
        memcpy(&x_lo, ladder->x + i, sizeof(x_lo));
        memcpy(&x_hi, ladder->x + i + DVEC_LEN, sizeof(x_hi));
        VEC phase = aov_dvec_pair_to_vecf(aov_fold_phase_dvec((DVEC){.data = x_lo.data * scale}), aov_fold_phase_dvec((DVEC){.data = x_hi.data * scale}));
        VEC weight = load_vecf_u(amp + i);
        VEC cos_phase = cos2pif_vec(phase);
        VEC sin_phase = sin2pif_vec(phase);
        store_vecf_u(state_r + i, (VEC){.data = weight.data * cos_phase.data});
        store_vecf_u(state_i + i, (VEC){.data = weight.data * sin_phase.data});
    }
    for (; i < n; ++i) {
        double phase = ladder->x[i] * scale;
        phase -= (double)(long long)phase;
        float ph = (float)phase;
        float weight = amp[i];
        state_r[i] = weight * cos2pif_tls(ph);
        state_i[i] = weight * sin2pif_tls(ph);
    }
}

/* Regenerate every harmonic at `block`. */
static inline void ladder_regenerate(twiddle_ladder *ladder, size_t block) {
    for (uint32_t h = 0; h < ladder->harmonics; ++h) ladder_regenerate_harmonic(ladder, h, block);
}

/* Level 0: state *= step for one harmonic (pad is a VEC_LEN multiple). */
static inline void ladder_apply_step(float *restrict state_r, float *restrict state_i, const float *restrict step_r, const float *restrict step_i, size_t pad) {
    for (size_t i = 0; i < pad; i += VEC_LEN) {
        VEC wr = {.data = *(const vecf_data *)(state_r + i)};
        VEC wi = {.data = *(const vecf_data *)(state_i + i)};
        VEC dr = {.data = *(const vecf_data *)(step_r + i)};
        VEC di = {.data = *(const vecf_data *)(step_i + i)};
        *(vecf_data *)(state_r + i) = wr.data * dr.data - wi.data * di.data;
        *(vecf_data *)(state_i + i) = wr.data * di.data + wi.data * dr.data;
    }
}

/* Advance the ladder to `block` (monotone during a sweep): resync from the exact
 * phase at an anchor, otherwise one complex multiply per harmonic. */
static inline void ladder_advance(twiddle_ladder *ladder, size_t block) {
    if (!ladder->step_r || (block % (size_t)ladder->anchor) == 0U) {
        ladder_regenerate(ladder, block);
        return;
    }
    const size_t pad = ladder->pad;
    for (uint32_t h = 0; h < ladder->harmonics; ++h) {
        const size_t offset = (size_t)h * pad;
        ladder_apply_step(ladder->state_r + offset, ladder->state_i + offset, ladder->step_r + offset, ladder->step_i + offset, pad);
    }
}

/* Pointer to one harmonic's current state (used as the NuFFT input directly). */
static inline void ladder_state(const twiddle_ladder *ladder, uint32_t harmonic, const float **real, const float **imag) {
    const size_t offset = (size_t)harmonic * ladder->pad;
    *real = ladder->state_r + offset;
    *imag = ladder->state_i + offset;
}

/* Build the per-harmonic read-only step (one complex array per harmonic) and the
 * block-0 state.  `dy`/`eps` are only used when `mult` is NULL (raw weight);
 * `y`/`ymean` are not needed here.  `num_blocks` decides whether the step is
 * needed at all. */
static inline int ladder_start(twiddle_ladder *ladder, const double *x, size_t n, size_t pad, uint32_t harmonics, uint32_t block_len, size_t num_blocks,
                               double fmin, double df, const float *mult, const float *dy, float eps) {
    if (!ladder || !x || n == 0U || pad < n || (pad % (size_t)VEC_LEN) != 0U) return -1;
    if (!double_is_finite_bits(fmin) || !double_is_finite_bits(df) || harmonics == 0U || block_len == 0U) return -1;

    const bool need_steps = num_blocks > 1U;
    if (ladder_alloc(ladder, harmonics, pad, need_steps) != 0) return -1;

    ladder->harmonics = harmonics;
    ladder->anchor = (uint32_t)MAX_TWIDDLE_REUSE;
    ladder->block_len = block_len;
    ladder->n = n;
    ladder->pad = pad;
    ladder->fmin = fmin;
    ladder->df = df;
    ladder->x = x;
    ladder->mult = mult;

    if (mult == NULL) {
        size_t i = 0;
        for (; i < n; ++i) ladder->weight[i] = 1.0f / ((dy[i] * dy[i]) + eps);
        for (; i < pad; ++i) ladder->weight[i] = 0.0f;
    }

    /* Zero state tail once: it is only ever multiplied afterwards. */
    for (uint32_t h = 0; h < harmonics; ++h) {
        const size_t offset = (size_t)h * pad;
        for (size_t i = n; i < pad; ++i) {
            ladder->state_r[offset + i] = 0.0f;
            ladder->state_i[offset + i] = 0.0f;
        }
    }

    if (need_steps) {
        for (uint32_t h = 0; h < harmonics; ++h) {
            const double scale = (double)(h + 1U) * df * (double)block_len;
            float *restrict step_r = ladder->step_r + (size_t)h * pad;
            float *restrict step_i = ladder->step_i + (size_t)h * pad;
            size_t i = 0;
            for (; i < n; ++i) {
                double phase = x[i] * scale;
                phase -= (double)(long long)phase;
                float ph = (float)phase;
                step_r[i] = cos2pif_tls(ph);
                step_i[i] = sin2pif_tls(ph);
            }
            for (; i < pad; ++i) { /* zero state tail * 1 stays zero */
                step_r[i] = 1.0f;
                step_i[i] = 0.0f;
            }
        }
    }

    ladder_regenerate(ladder, 0);
    return 0;
}

/* Trig sums of one block: unweighted (q = 1..Q) and y-weighted (q = 1..Q/2)
 * harmonic sums, written with row stride `count`.  Both share one ladder state
 * and one advance per (block, q); the weighted input is (y - ymean) * state.
 * `tmp_r`/`tmp_i` are block-sized scratch used when the rows are not exactly
 * block_len long (partial last block).
 *
 * `slot_order` asks for the deferred bit reversal: the NuFFT hands back the
 * block in slot order (kept FFT positions only), the rows are block_len long
 * either way, and the caller turns the *results* into grid order afterwards
 * with nufft1_results_to_grid()/nufft1_slots_to_grid().  Full blocks only -
 * a partial last block must stay frequency-ordered. */
static inline void ladder_block_sums(twiddle_ladder *ladder, nufft1_workspace *workspace, const float *y, float ymean, float *Sw, float *Cw, float *Syw,
                                     float *Cyw, float *tmp_r, float *tmp_i, uint32_t count, bool slot_order) {
    const size_t n = ladder->n;
    const size_t pad = ladder->pad;
    const size_t n_vec = n - (n % (size_t)VEC_LEN);
    const size_t stride = slot_order ? (size_t)ladder->block_len : (size_t)count;
    const uint32_t weighted_max = ladder->harmonics / 2U;
    const size_t bytes = stride * sizeof(float);
    const bool direct = slot_order || (stride == (size_t)ladder->block_len); /* row is exactly block_len long */
    const VEC ymean_vec = SET_VEC(ymean);

    for (uint32_t q = 1; q <= ladder->harmonics; ++q) {
        const size_t offset = (size_t)(q - 1U) * pad;
        const float *state_r = ladder->state_r + offset;
        const float *state_i = ladder->state_i + offset;
        /* Full blocks write the NuFFT output straight into the sum rows. */
        float *sum_c = direct ? Cw + (size_t)q * stride : tmp_r;
        float *sum_s = direct ? Sw + (size_t)q * stride : tmp_i;

        if (slot_order)
            nufft1_execute_slots(workspace, state_r, state_i, sum_c, sum_s, (int)q);
        else
            nufft1_execute(workspace, state_r, state_i, sum_c, sum_s, (int)q);
        if (!direct) {
            memcpy(Cw + (size_t)q * stride, tmp_r, bytes);
            memcpy(Sw + (size_t)q * stride, tmp_i, bytes);
        }

        if (q <= weighted_max) {
            size_t i = 0;
            for (; i < n_vec; i += VEC_LEN) {
                VEC yv = load_vecf_u(y + i);
                VEC sr = load_vecf_u(state_r + i);
                VEC si = load_vecf_u(state_i + i);
                VEC scale = {.data = yv.data - ymean_vec.data};
                store_vecf_u(ladder->yin_r + i, (VEC){.data = scale.data * sr.data});
                store_vecf_u(ladder->yin_i + i, (VEC){.data = scale.data * si.data});
            }
            for (; i < n; ++i) {
                float scale = y[i] - ymean;
                ladder->yin_r[i] = scale * state_r[i];
                ladder->yin_i[i] = scale * state_i[i];
            }
            float *wy_c = direct ? Cyw + (size_t)q * stride : tmp_r;
            float *wy_s = direct ? Syw + (size_t)q * stride : tmp_i;
            if (slot_order)
                nufft1_execute_slots(workspace, ladder->yin_r, ladder->yin_i, wy_c, wy_s, (int)q);
            else
                nufft1_execute(workspace, ladder->yin_r, ladder->yin_i, wy_c, wy_s, (int)q);
            if (!direct) {
                memcpy(Cyw + (size_t)q * stride, tmp_r, bytes);
                memcpy(Syw + (size_t)q * stride, tmp_i, bytes);
            }
        }
    }
}

#endif  // LADDER_H
