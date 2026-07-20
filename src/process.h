#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "utils/readout.h"

#include "profile.h"
#include "utils/simd.h"
#include "utils/strings.h"
#include "utils/trig.h"

static inline float correctPowerTmp(float K, float nInv) {
    VEC K_tmp = SET_VEC(0.0f);
    K_tmp.values[0] = K;
    K_tmp = correctPower(K_tmp, nInv);
    return K_tmp.values[0];
}

static inline double effective_grid_fmin(const parameters *params, double delta_t) {
    if (params->fmin == 0.0f && delta_t > 0.0) return 2.0 / delta_t;
    return (double)params->fmin;
}

static inline uint32_t target_frequency_count(const parameters *params, double fmin, double fstep) {
    if (fstep <= 0.0 || (double)params->fmax <= fmin) return 0;
    return (uint32_t)floor(((double)params->fmax - fmin) / fstep) + 1U;
}

static inline bool activate_target_nufft_plan(buffer_t *buffer, const parameters *params, uint32_t nfreq) {
    if (!buffer || !params || params->nufftPlanCount == 0U || !buffer->nufftWorkspace) return false;

    uint32_t plan_idx = params->nufftPlanCount - 1U;
    int target_nfreq = nfreq > 0U ? (int)nfreq : 1;
    if ((uint64_t)buffer->n * 2U < (uint64_t)params->maxLen || (uint64_t)nfreq * 2U < (uint64_t)params->maxFreqCount) {
        double beta = params->gridMode == NUFFT1_PSWF21 ? F_PSWF21_BETA : F_PSWF43_BETA;
        double gamma = params->gridMode == NUFFT1_PSWF21 ? F_PSWF21_GAMMA : F_PSWF43_GAMMA;
        int plan_degree = periodogram_uses_aov(params->periodogramMethod) ? params->nterms : 0;
        int base_len = nufft1_optimize_plan_size(target_nfreq, (int)buffer->n, plan_degree, F_ALPHA, beta, gamma, pswf_backend(params->gridMode));
        if (params->gridMode == NUFFT1_PSWF43) base_len = nufft1_pswf43_plan_len_from_base(base_len);
        if (base_len < (1 << 8)) base_len = 1 << 8;
        if ((uint32_t)base_len > params->gridLen) base_len = (int)params->gridLen;
        if (params->gridMode == NUFFT1_PSWF43) base_len = nufft1_pswf43_plan_len_from_base(base_len);
        if ((uint32_t)base_len > params->gridLen) base_len = (int)params->gridLen;

        for (uint32_t i = 0; i < params->nufftPlanCount; ++i) {
            if (params->nufftPlanCache[i].gridLen >= (uint32_t)base_len) {
                plan_idx = i;
                break;
            }
        }
    }

    const nufft_plan_cache_entry_t *entry = &params->nufftPlanCache[plan_idx];
    if (nufft1_workspace_set_plan(buffer->nufftWorkspace, entry->plan) != NUFFT1_UTIL_OK) return false;

    buffer->activePlanIndex = plan_idx;
    buffer->activeGridLen = entry->gridLen;
    buffer->activeOutputLen = entry->outputLen;
    buffer->activeLadderLevels = (uint32_t)nufft1_twiddle_ladder_levels(target_nfreq, (int)entry->outputLen);
    if (buffer->activeLadderLevels > NUFFT_LADDER_LEVEL_CAP) buffer->activeLadderLevels = NUFFT_LADDER_LEVEL_CAP;
    if (buffer->activeLadderLevels == 0U) buffer->activeLadderLevels = 1U;
    return true;
}

static inline void fill_complex_input(buffer_t *buffer, double base_freq, double freq_factor) {
    uint32_t i = 0;
    for (; i < buffer->n; ++i) {
        float phase = (float)(buffer->x[i] * base_freq * freq_factor);
        buffer->inputReal[i] = buffer->wy[i] * cos2pif_tls(phase);
        buffer->inputImag[i] = buffer->wy[i] * sin2pif_tls(phase);
    }
    for (; i < buffer->paddedLen; ++i) {
        buffer->inputReal[i] = 0.0f;
        buffer->inputImag[i] = 0.0f;
    }
}

static inline void fill_twiddle_ladder(buffer_t *buffer, double fstep, double freq_factor) {
    for (uint32_t level = 0; level < buffer->activeLadderLevels; ++level) {
        size_t offset = (size_t)level * (size_t)buffer->paddedLen;
        double advance = nufft1_twiddle_ladder_advance((int)buffer->activeOutputLen, (int)level);
        uint32_t i = 0;
        for (; i < buffer->n; ++i) {
            float phase = (float)(buffer->x[i] * fstep * freq_factor * advance);
            buffer->deltaReal[offset + i] = cos2pif_tls(phase);
            buffer->deltaImag[offset + i] = sin2pif_tls(phase);
        }
        for (; i < buffer->paddedLen; ++i) {
            buffer->deltaReal[offset + i] = 1.0f;
            buffer->deltaImag[offset + i] = 0.0f;
        }
    }
}

static inline void fill_twiddle_ladder_backward(buffer_t *buffer, double fstep, double freq_factor) {
    for (uint32_t level = 0; level < buffer->activeLadderLevels; ++level) {
        size_t offset = (size_t)level * (size_t)buffer->paddedLen;
        double advance = nufft1_twiddle_ladder_advance((int)buffer->activeOutputLen, (int)level);
        uint32_t i = 0;
        for (; i < buffer->n; ++i) {
            float phase = (float)(buffer->x[i] * fstep * freq_factor * (-advance));
            buffer->deltaBackReal[offset + i] = cos2pif_tls(phase);
            buffer->deltaBackImag[offset + i] = sin2pif_tls(phase);
        }
        for (; i < buffer->paddedLen; ++i) {
            buffer->deltaBackReal[offset + i] = 1.0f;
            buffer->deltaBackImag[offset + i] = 0.0f;
        }
    }
}

// Compute work state at arbitrary block K from pristine input using hierarchical ladder decomposition
// O(log_8(K) * 7) complex multiplications - bounded error
// Each ladder level L is set to the state at the last carry point for that level:
//   work[L] = input * prod_{j>=L} delta[j]^d[j]  where d[j] are base-8 digits of block_idx
static inline void compute_work_at_block(buffer_t *buffer, size_t block_idx) {
    // Start with pristine input at level 0
    memcpy(buffer->workReal, buffer->inputReal, (size_t)buffer->paddedLen * sizeof(float));
    memcpy(buffer->workImag, buffer->inputImag, (size_t)buffer->paddedLen * sizeof(float));

    if (block_idx == 0) {
        // Copy to all levels
        for (uint32_t level = 1; level < buffer->activeLadderLevels; ++level) {
            size_t offset = (size_t)level * (size_t)buffer->paddedLen;
            memcpy(buffer->workReal + offset, buffer->inputReal, (size_t)buffer->paddedLen * sizeof(float));
            memcpy(buffer->workImag + offset, buffer->inputImag, (size_t)buffer->paddedLen * sizeof(float));
        }
        return;
    }

    // Extract base-8 digits
    size_t digits[NUFFT_LADDER_LEVEL_CAP];
    size_t remaining = block_idx;
    uint32_t max_level = 0;
    for (uint32_t level = 0; level < buffer->activeLadderLevels; ++level) {
        digits[level] = remaining % 8U;
        remaining /= 8U;
        if (digits[level] > 0) max_level = level;
    }

    // Process levels from highest to lowest.
    // At each level L, multiply work[0] by delta[L]^digits[L],
    // then save the result to work[L] (the carry-point state for this level).
    // work[L] = input * prod_{j>=L} delta[j]^d[j]
    uint32_t i, n_vec;
    for (int level = (int)max_level; level >= 0; --level) {
        // Multiply work[0] by delta[level] exactly digits[level] times
        if (digits[level] > 0) {
            size_t offset = (size_t)level * (size_t)buffer->paddedLen;
            for (size_t d = 0; d < digits[level]; ++d) {
                n_vec = buffer->paddedLen - (buffer->paddedLen % VEC_LEN);
                for (i = 0; i < n_vec; i += VEC_LEN) {
                    VEC wr = {.data = *(vecf_data *)(buffer->workReal + i)};
                    VEC wi = {.data = *(vecf_data *)(buffer->workImag + i)};
                    VEC dr = {.data = *(vecf_data *)(buffer->deltaReal + offset + i)};
                    VEC di = {.data = *(vecf_data *)(buffer->deltaImag + offset + i)};
                    *(vecf_data *)(buffer->workReal + i) = wr.data * dr.data - wi.data * di.data;
                    *(vecf_data *)(buffer->workImag + i) = wr.data * di.data + wi.data * dr.data;
                }
                for (; i < buffer->paddedLen; ++i) {
                    float wr = buffer->workReal[i];
                    float wi = buffer->workImag[i];
                    float dr = buffer->deltaReal[offset + i];
                    float di = buffer->deltaImag[offset + i];
                    buffer->workReal[i] = wr * dr - wi * di;
                    buffer->workImag[i] = wr * di + wi * dr;
                }
            }
        }

        // Save current state to work[level] AFTER multiplication (carry-point state)
        if (level > 0) {
            size_t dst_offset = (size_t)level * (size_t)buffer->paddedLen;
            memcpy(buffer->workReal + dst_offset, buffer->workReal, (size_t)buffer->paddedLen * sizeof(float));
            memcpy(buffer->workImag + dst_offset, buffer->workImag, (size_t)buffer->paddedLen * sizeof(float));
        }
    }

    // Fill any remaining higher levels with input (they have digit=0)
    for (uint32_t level = max_level + 1; level < buffer->activeLadderLevels; ++level) {
        size_t offset = (size_t)level * (size_t)buffer->paddedLen;
        memcpy(buffer->workReal + offset, buffer->inputReal, (size_t)buffer->paddedLen * sizeof(float));
        memcpy(buffer->workImag + offset, buffer->inputImag, (size_t)buffer->paddedLen * sizeof(float));
    }
}

static inline void reset_work_ladder(buffer_t *buffer) {
    for (uint32_t level = 0; level < buffer->activeLadderLevels; ++level) {
        size_t offset = (size_t)level * (size_t)buffer->paddedLen;
        memcpy(buffer->workReal + offset, buffer->inputReal, (size_t)buffer->paddedLen * sizeof(float));
        memcpy(buffer->workImag + offset, buffer->inputImag, (size_t)buffer->paddedLen * sizeof(float));
    }
}

static inline void advance_work_ladder(buffer_t *buffer, size_t next_block) {
    int level = nufft1_twiddle_ladder_carry_level(next_block, (int)buffer->activeLadderLevels);
    if (level >= (int)buffer->activeLadderLevels) level = (int)buffer->activeLadderLevels - 1;
    size_t offset = (size_t)level * (size_t)buffer->paddedLen;

    uint32_t i = 0;
    uint32_t n_vec = buffer->paddedLen - (buffer->paddedLen % VEC_LEN);
    for (; i < n_vec; i += VEC_LEN) {
        VEC wr = {.data = *(vecf_data *)(buffer->workReal + offset + i)};
        VEC wi = {.data = *(vecf_data *)(buffer->workImag + offset + i)};
        VEC dr = {.data = *(vecf_data *)(buffer->deltaReal + offset + i)};
        VEC di = {.data = *(vecf_data *)(buffer->deltaImag + offset + i)};
        *(vecf_data *)(buffer->workReal + offset + i) = wr.data * dr.data - wi.data * di.data;
        *(vecf_data *)(buffer->workImag + offset + i) = wr.data * di.data + wi.data * dr.data;
    }
    for (; i < buffer->paddedLen; ++i) {
        float wr = buffer->workReal[offset + i];
        float wi = buffer->workImag[offset + i];
        float dr = buffer->deltaReal[offset + i];
        float di = buffer->deltaImag[offset + i];
        buffer->workReal[offset + i] = wr * dr - wi * di;
        buffer->workImag[offset + i] = wr * di + wi * dr;
    }

    for (int dst_level = level; dst_level > 0; --dst_level) {
        size_t dst = (size_t)(dst_level - 1) * (size_t)buffer->paddedLen;
        size_t src = (size_t)dst_level * (size_t)buffer->paddedLen;
        memcpy(buffer->workReal + dst, buffer->workReal + src, (size_t)buffer->paddedLen * sizeof(float));
        memcpy(buffer->workImag + dst, buffer->workImag + src, (size_t)buffer->paddedLen * sizeof(float));
    }
}

static inline void accumulate_power(buffer_t *buffer, uint32_t base, uint32_t count) {
    uint32_t i = 0;
    uint32_t n_vec = count - (count % VEC_LEN);
    for (; i < n_vec; i += VEC_LEN) {
        VEC r = {.data = *(vecf_data *)(buffer->blockReal + i)};
        VEC im = {.data = *(vecf_data *)(buffer->blockImag + i)};
        VEC p = {.data = *(vecf_data *)(buffer->power + base + i)};
        *(vecf_data *)(buffer->power + base + i) = p.data + r.data * r.data + im.data * im.data;
    }
    for (; i < count; ++i) {
        buffer->power[base + i] += (buffer->blockReal[i] * buffer->blockReal[i]) + (buffer->blockImag[i] * buffer->blockImag[i]);
    }
}

static inline void execute_nufft_sweep(buffer_t *buffer, parameters *params, double fmin, double fstep, uint32_t nfreq) {
    (void)params;
    memset(buffer->power, 0, ((size_t)nfreq + 1U) * sizeof(float));
    {
        PROFILE_START(precompute);
        nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);
        PROFILE_END(precompute, PROFILE_PHASE_PRECOMPUTE);
    }

    size_t num_blocks = ((size_t)nfreq + (size_t)buffer->activeOutputLen - 1U) / (size_t)buffer->activeOutputLen;
    for (int t = 0; t < buffer->terms; ++t) {
        double freq_factor = (double)(t + 1);
        {
            PROFILE_START(fill_input);
            fill_complex_input(buffer, fmin, freq_factor);
            PROFILE_END(fill_input, PROFILE_PHASE_FILL_INPUT);
        }
        {
            PROFILE_START(ladder_setup);
            fill_twiddle_ladder(buffer, fstep, freq_factor);
            PROFILE_END(ladder_setup, PROFILE_PHASE_LADDER_SETUP);
        }
        {
            PROFILE_START(reset_ladder);
            reset_work_ladder(buffer);
            PROFILE_END(reset_ladder, PROFILE_PHASE_RESET_LADDER);
        }

        for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
            uint32_t base = (uint32_t)(block_idx * (size_t)buffer->activeOutputLen);
            uint32_t count = buffer->activeOutputLen;
            if (base + count > nfreq) count = nfreq - base;
            {
                PROFILE_START(nufft_execute);
                nufft1_execute(buffer->nufftWorkspace, buffer->workReal, buffer->workImag, buffer->blockReal, buffer->blockImag, t + 1);
                PROFILE_END(nufft_execute, PROFILE_PHASE_NUFFT_EXECUTE);
            }
            {
                PROFILE_START(power_accumulate);
                accumulate_power(buffer, base, count);
                PROFILE_END(power_accumulate, PROFILE_PHASE_POWER_ACCUMULATE);
            }
            if (block_idx + 1U < num_blocks) {
                PROFILE_START(ladder_advance);
                advance_work_ladder(buffer, block_idx + 1U);
                PROFILE_END(ladder_advance, PROFILE_PHASE_LADDER_ADVANCE);
            }
        }
    }
}

// ======================== Frequency Range Slicing (Parallel Sweep) ========================
// Divides [fmin, fmax] into N contiguous sub-ranges (one per thread).
// Each worker runs the unmodified serial sweep on its sub-range, then results are merged.
// Adjacent sub-ranges overlap by exactly 1 frequency for seamless peak detection.
// Minor FP differences vs full serial (lower accumulated ladder error per worker).

typedef struct {
    uint32_t freq_start;   // global start index
    uint32_t slice_nfreq;  // frequencies in this slice
    int buffer_idx;        // index into params->buffers (set by worker; may differ from workset index due to work stealing)
    int status;
} freq_slice_workset_t;

typedef struct {
    buffer_t *primary;
    parameters *params;
    double fmin;
    double fstep;
    uint32_t nfreq;
    freq_slice_workset_t *worksets;
} freq_slice_context_t;

static inline void freq_slice_ihs_worker(void *data, long i, int thread_id) {
    freq_slice_context_t *ctx = (freq_slice_context_t *)data;
    freq_slice_workset_t *work = &ctx->worksets[i];
    buffer_t *buf = ctx->params->buffers[thread_id];
    buffer_t *primary = ctx->primary;

    work->buffer_idx = thread_id;

    if (work->slice_nfreq == 0) {
        work->status = 0;
        return;
    }

    // Save worker's own allocations before sharing primary's data pointers
    double *own_x = buf->x;
    float *own_y = buf->y;
    float *own_dy = buf->dy;
    float *own_wy = buf->wy;

    // Sync worker buffer metadata from primary
    buf->n = primary->n;
    buf->paddedLen = primary->paddedLen;
    buf->terms = primary->terms;
    buf->x = primary->x;
    buf->y = primary->y;
    buf->dy = primary->dy;
    buf->wy = primary->wy;

    // Use the same plan as primary (ensures identical block size and output length)
    const nufft_plan_cache_entry_t *entry = &ctx->params->nufftPlanCache[primary->activePlanIndex];
    if (nufft1_workspace_set_plan(buf->nufftWorkspace, entry->plan) != NUFFT1_UTIL_OK) {
        work->status = -1;
        buf->x = own_x;
        buf->y = own_y;
        buf->dy = own_dy;
        buf->wy = own_wy;
        return;
    }
    buf->activePlanIndex = primary->activePlanIndex;
    buf->activeGridLen = primary->activeGridLen;
    buf->activeOutputLen = primary->activeOutputLen;
    buf->activeLadderLevels = primary->activeLadderLevels;

    // Run serial IHS sweep on this worker's frequency sub-range
    double slice_fmin = ctx->fmin + (double)work->freq_start * ctx->fstep;
    execute_nufft_sweep(buf, ctx->params, slice_fmin, ctx->fstep, work->slice_nfreq);
    work->status = 0;

    // Restore worker's own pointers so free_buffer doesn't free primary's data
    buf->x = own_x;
    buf->y = own_y;
    buf->dy = own_dy;
    buf->wy = own_wy;
}

static inline bool execute_nufft_sweep_parallel(buffer_t *buffer, parameters *params, kt_forpool_t *pool, double fmin, double fstep, uint32_t nfreq) {
    if (nfreq == 0) return true;
    if (!pool || pool->n_threads < 2 || !params->buffers || params->nbuffers < pool->n_threads) {
        execute_nufft_sweep(buffer, params, fmin, fstep, nfreq);
        return true;
    }

    int nthreads = pool->n_threads;
    uint32_t nworksets = (uint32_t)nthreads;
    if (nworksets > nfreq) nworksets = nfreq;
    if (nworksets < 2U) nworksets = 2U;

    // Frequency range slicing with 1-frequency overlap at boundaries
    uint32_t freqs_per_worker = (nfreq + nworksets - 2U) / (nworksets - 1U);

    freq_slice_workset_t *worksets = calloc(nworksets, sizeof(*worksets));
    if (!worksets) return false;

    for (uint32_t w = 0; w < nworksets; ++w) {
        uint32_t start = w * (freqs_per_worker - 1U);
        uint32_t count = freqs_per_worker;
        if (start >= nfreq) {
            start = nfreq;
            count = 0;
        } else if (start + count > nfreq) {
            count = nfreq - start;
        }
        worksets[w].freq_start = start;
        worksets[w].slice_nfreq = count;
        worksets[w].status = -1;
    }

    freq_slice_context_t ctx = {.primary = buffer, .params = params, .fmin = fmin, .fstep = fstep, .nfreq = nfreq, .worksets = worksets};

    kt_forpool(pool, freq_slice_ihs_worker, &ctx, (long)nworksets);

    bool ok = true;
    for (uint32_t w = 0; w < nworksets; ++w) {
        if (worksets[w].status != 0) ok = false;
    }

    // Merge: copy each worker's power into primary, skipping overlap for workers > 0
    if (ok) {
        for (uint32_t w = 0; w < nworksets; ++w) {
            if (worksets[w].slice_nfreq == 0) continue;
            buffer_t *buf = params->buffers[worksets[w].buffer_idx];
            uint32_t skip = (w == 0) ? 0U : 1U;
            memcpy(buffer->power + worksets[w].freq_start + skip, buf->power + skip, (size_t)(worksets[w].slice_nfreq - skip) * sizeof(float));
        }
    }

    free(worksets);
    return ok;
}

#include "utils/aov.h"

// ======================== Parallel AoV Sweep (Frequency Range Slicing) ========================

static inline void freq_slice_aov_worker(void *data, long i, int thread_id) {
    freq_slice_context_t *ctx = (freq_slice_context_t *)data;
    freq_slice_workset_t *work = &ctx->worksets[i];
    buffer_t *buf = ctx->params->buffers[thread_id];
    buffer_t *primary = ctx->primary;
    parameters *params = ctx->params;

    work->buffer_idx = thread_id;

    if (work->slice_nfreq == 0) {
        work->status = 0;
        return;
    }

    // Save worker's own allocations before sharing primary's data pointers
    double *own_x = buf->x;
    float *own_y = buf->y;
    float *own_dy = buf->dy;
    float *own_wy = buf->wy;

    // Sync worker buffer metadata from primary
    buf->n = primary->n;
    buf->paddedLen = primary->paddedLen;
    buf->terms = primary->terms;
    buf->neff = primary->neff;
    buf->amp_neff = primary->amp_neff;
    buf->x = primary->x;
    buf->y = primary->y;
    buf->dy = primary->dy;
    buf->wy = primary->wy;

    // Use the same plan as primary (ensures identical block size and output length)
    const nufft_plan_cache_entry_t *entry = &params->nufftPlanCache[primary->activePlanIndex];
    if (nufft1_workspace_set_plan(buf->nufftWorkspace, entry->plan) != NUFFT1_UTIL_OK) {
        work->status = -1;
        buf->x = own_x;
        buf->y = own_y;
        buf->dy = own_dy;
        buf->wy = own_wy;
        return;
    }
    buf->activePlanIndex = primary->activePlanIndex;
    buf->activeGridLen = primary->activeGridLen;
    buf->activeOutputLen = primary->activeOutputLen;
    buf->activeLadderLevels = primary->activeLadderLevels;

    // Run serial AoV sweep on this worker's frequency sub-range (store_power=true, no spectrum/peaks)
    double slice_fmin = ctx->fmin + (double)work->freq_start * ctx->fstep;
    char stringBuff[64];
    int status = execute_aov_sweep(buf, params, slice_fmin, ctx->fstep, work->slice_nfreq, INFINITY, 0.0, false, true, false, 0, stringBuff);
    work->status = (status == 0) ? 0 : -1;

    // Restore worker's own pointers so free_buffer doesn't free primary's data
    buf->x = own_x;
    buf->y = own_y;
    buf->dy = own_dy;
    buf->wy = own_wy;
}

static inline bool execute_aov_sweep_parallel(buffer_t *buffer, parameters *params, kt_forpool_t *pool, double fmin, double fstep, uint32_t nfreq) {
    if (nfreq == 0) return true;
    if (!pool || pool->n_threads < 2 || !params->buffers || params->nbuffers < pool->n_threads) {
        // Fallback to serial
        char stringBuff[64];
        return execute_aov_sweep(buffer, params, fmin, fstep, nfreq, INFINITY, 0.0, false, true, false, 0, stringBuff) == 0;
    }

    // Handle edge cases identically to serial version
    if (!aov_target_has_dof(buffer, params->nterms)) {
        memset(buffer->power, 0, (size_t)nfreq * sizeof(float));
        return true;
    }
    aov_reference_t ref;
    if (!aov_prepare_reference(buffer, params->epsilon, &ref)) {
        memset(buffer->power, 0, (size_t)nfreq * sizeof(float));
        return true;
    }

    int nthreads = pool->n_threads;
    uint32_t nworksets = (uint32_t)nthreads;
    if (nworksets > nfreq) nworksets = nfreq;
    if (nworksets < 2U) nworksets = 2U;

    // Frequency range slicing with 1-frequency overlap at boundaries
    uint32_t freqs_per_worker = (nfreq + nworksets - 2U) / (nworksets - 1U);

    freq_slice_workset_t *worksets = calloc(nworksets, sizeof(*worksets));
    if (!worksets) return false;

    for (uint32_t w = 0; w < nworksets; ++w) {
        uint32_t start = w * (freqs_per_worker - 1U);
        uint32_t count = freqs_per_worker;
        if (start >= nfreq) {
            start = nfreq;
            count = 0;
        } else if (start + count > nfreq) {
            count = nfreq - start;
        }
        worksets[w].freq_start = start;
        worksets[w].slice_nfreq = count;
        worksets[w].status = -1;
    }

    freq_slice_context_t ctx = {.primary = buffer, .params = params, .fmin = fmin, .fstep = fstep, .nfreq = nfreq, .worksets = worksets};

    kt_forpool(pool, freq_slice_aov_worker, &ctx, (long)nworksets);

    bool ok = true;
    for (uint32_t w = 0; w < nworksets; ++w) {
        if (worksets[w].status != 0) ok = false;
    }

    // Merge: copy each worker's power into primary, skipping overlap for workers > 0
    if (ok) {
        for (uint32_t w = 0; w < nworksets; ++w) {
            if (worksets[w].slice_nfreq == 0) continue;
            buffer_t *buf = params->buffers[worksets[w].buffer_idx];
            uint32_t skip = (w == 0) ? 0U : 1U;
            memcpy(buffer->power + worksets[w].freq_start + skip, buf->power + skip, (size_t)(worksets[w].slice_nfreq - skip) * sizeof(float));
        }
    }

    free(worksets);
    return ok;
}

static inline bool quadratic_peak_position(double freq, double fstep, float left, float center, float right, float oversampling, double *peak_freq,
                                           float *peak_magnitude) {
    if (!float_is_finite_bits(left) || !float_is_finite_bits(center) || !float_is_finite_bits(right)) return false;
    *peak_freq = freq;
    *peak_magnitude = center;
    if (oversampling < 3.0f) return true;

    double x0 = freq - fstep;
    double x1 = freq;
    double x2 = freq + fstep;
    double slope01 = ((double)center - (double)left) / (x1 - x0);
    double slope12 = ((double)right - (double)center) / (x2 - x1);
    double curvature = (slope12 - slope01) / (x2 - x0);
    if (curvature == 0.0) return true;

    double linear = slope01 - curvature * (x0 + x1);
    double vx = -linear / (2.0 * curvature);
    double vy = (double)left + slope01 * (vx - x0) + curvature * (vx - x0) * (vx - x1);
    float vyf = (float)vy;
    if (double_is_finite_bits(vx) && float_is_finite_bits(vyf) && vx >= x0 && vx <= x2) {
        *peak_freq = vx;
        *peak_magnitude = vyf;
    }
    return true;
}

static inline void capture_spectrum_column(buffer_t *buffer, const parameters *params, const eval_method_t *eval_method, bool use_aov, bool direct_eval_grid,
                                           int aov_n_eff, double fmin, double fstep, uint32_t nfreq, int nterms, float *spectrum_matrix, uint32_t column) {
    if (!spectrum_matrix) return;
    float *dst = spectrum_matrix + ((size_t)column * (size_t)nfreq);
    for (uint32_t i = 0; i < nfreq; ++i) {
        double freq = fmin + ((double)i * fstep);
        float magnitude;
        if (mode_uses_direct_eval_grid(params->mode)) {
            if (direct_eval_grid) {
                magnitude = buffer->power ? buffer->power[i] : get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL);
            } else if (use_aov) {
                float r2 = buffer->power ? buffer->power[i] : aov_get_stat(buffer, params, freq);
                magnitude = aov_likelihood_from_r2(r2, nterms, aov_n_eff);
            } else {
                magnitude = buffer->power ? buffer->power[i] : get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL);
            }
        } else if (use_aov) {
            float r2 = buffer->power[i];
            magnitude = aov_likelihood_from_r2(r2, nterms, aov_n_eff);
        } else {
            magnitude = correct_ihs_res(buffer->power[i], nterms);
        }
        dst[i] = float_is_finite_bits(magnitude) ? magnitude : 0.0f;
    }
}

typedef struct {
    uint32_t core_begin;
    uint32_t core_end;
    uint32_t eval_begin;
    uint32_t eval_end;
    float *values;
    peak_t *peaks;
    uint32_t nPeaks;
    sds spectrum;
    int status;
} direct_grid_workset_t;

typedef struct {
    buffer_t *buffer;
    parameters *params;
    const eval_method_t *eval_method;
    bool threshold_valid;
    bool scan_peaks;
    bool write_spectrum;
    float threshold;
    double fmin;
    double fstep;
    double df;
    uint32_t nfreq;
    int precision;
    buffer_t *scratch;
    direct_grid_workset_t *worksets;
} direct_grid_context_t;

static inline uint32_t direct_grid_chunk_size(uint32_t n) {
    uint32_t chunk = 1U;
    uint64_t product = (uint64_t)n;
    while (product <= (UINT64_C(1) << 20) && chunk < (UINT32_MAX >> 1U)) {
        chunk <<= 1U;
        product = (uint64_t)chunk * (uint64_t)n;
    }
    if (chunk < 128U) chunk = 128U;
    return chunk;
}

static inline float direct_grid_scan_value(buffer_t *buffer, const parameters *params, const eval_method_t *eval_method, double freq) {
    return get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL);
}

static inline float direct_grid_spectrum_value(float scan_value) { return float_is_finite_bits(scan_value) ? scan_value : 0.0f; }

static inline int alloc_direct_scratch(buffer_t *scratch, const buffer_t *source, const parameters *params, bool use_aov) {
    memset(scratch, 0, sizeof(*scratch));
    scratch->allocated = true;
    scratch->len = source->len;
    scratch->n = source->n;
    scratch->terms = source->terms;
    scratch->neff = source->neff;
    scratch->amp_neff = source->amp_neff;
    scratch->magnitude = source->magnitude;
    scratch->maxFreqCount = source->maxFreqCount;
    scratch->x = source->x;
    scratch->y = source->y;
    scratch->dy = source->dy;
    scratch->wy = source->wy;
    scratch->pidx = (size_t *)malloc(1024U * sizeof(size_t));
    if (!scratch->pidx) return -1;

    scratch->buf = calloc(3U, sizeof(void *));
    if (!scratch->buf) return -1;
    for (int i = 0; i < 3; ++i) {
        scratch->buf[i] = aligned_alloc(64, round_buffer((size_t)params->maxLen * sizeof(uint64_t)));
        if (!scratch->buf[i]) return -1;
    }

    if (use_aov) {
        scratch->aovScratchLen = params->nterms > 0 ? (size_t)(26 * params->nterms + 16) : 0U;
        if (scratch->aovScratchLen > 0U) {
            scratch->aovScratch = aligned_alloc(64, round_buffer(scratch->aovScratchLen * sizeof(float)));
            if (!scratch->aovScratch) return -1;
        }
    }
    return 0;
}

static inline void free_direct_scratch(buffer_t *scratch) {
    free(scratch->pidx);
    scratch->pidx = NULL;
    if (scratch->buf) {
        for (int i = 0; i < 3; ++i) {
            free(scratch->buf[i]);
            scratch->buf[i] = NULL;
        }
        free(scratch->buf);
        scratch->buf = NULL;
    }
    free(scratch->aovScratch);
    scratch->aovScratch = NULL;
    scratch->aovScratchLen = 0U;
    scratch->x = NULL;
    scratch->y = NULL;
    scratch->dy = NULL;
    scratch->wy = NULL;
    scratch->peaks = NULL;
    scratch->nPeaks = 0U;
    scratch->allocated = false;
}

static inline float direct_merge_rank(const parameters *params, const eval_method_t *method, const peak_t *peak) {
    return mode_defers_peak_evaluation(params->mode) ? peak->p : eval_peak_rank(method, peak);
}

static inline void direct_merge_peak_stable(peak_t *dst, uint32_t *dst_count, int max_peaks, const peak_t *peak, const parameters *params,
                                            const eval_method_t *method) {
    if (max_peaks <= 0 || !peak) return;
    float rank = direct_merge_rank(params, method, peak);
    if (!float_is_finite_bits(rank)) return;
    uint32_t count = *dst_count;
    if (count >= (uint32_t)max_peaks) {
        float last_rank = direct_merge_rank(params, method, &dst[max_peaks - 1]);
        if (rank <= last_rank) return;
    }

    uint32_t idx = count;
    while (idx > 0U && rank > direct_merge_rank(params, method, &dst[idx - 1U])) {
        --idx;
    }
    if (count < (uint32_t)max_peaks) {
        ++count;
        *dst_count = count;
    }
    for (uint32_t i = count - 1U; i > idx; --i) dst[i] = dst[i - 1U];
    if (idx < (uint32_t)max_peaks) dst[idx] = *peak;
}

static inline void direct_grid_worker(void *data, long i, int thread_id) {
    direct_grid_context_t *ctx = (direct_grid_context_t *)data;
    direct_grid_workset_t *work = &ctx->worksets[i];
    buffer_t *scratch = &ctx->scratch[thread_id];
    char stringBuff[64];
    scratch->peaks = work->peaks;
    scratch->nPeaks = 0U;

    for (uint32_t idx = work->eval_begin; idx < work->eval_end; ++idx) {
        double freq = ctx->fmin + ((double)idx * ctx->fstep);
        float value = direct_grid_scan_value(scratch, ctx->params, ctx->eval_method, freq);
        work->values[idx - work->eval_begin] = value;
        if (idx >= work->core_begin && idx < work->core_end) {
            ctx->buffer->power[idx] = value;
            if (ctx->write_spectrum) {
                appendFreq(freq, direct_grid_spectrum_value(value), ctx->precision, ctx->params->outputPeriod, &work->spectrum, stringBuff);
            }
        }
    }

    if (ctx->scan_peaks) {
        for (uint32_t idx = work->core_begin; idx < work->core_end; ++idx) {
            if (idx == 0U || idx + 1U >= ctx->nfreq) continue;
            float left = work->values[(idx - 1U) - work->eval_begin];
            float center = work->values[idx - work->eval_begin];
            float right = work->values[(idx + 1U) - work->eval_begin];
            double freq = ctx->fmin + ((double)idx * ctx->fstep);
            double peak_freq = freq;
            float peak_magnitude = center;
            if (!quadratic_peak_position(freq, ctx->fstep, left, center, right, ctx->params->oversamplingFactor, &peak_freq, &peak_magnitude)) continue;
            if (center > left && center > right && (center > ctx->threshold || (ctx->threshold_valid && mode_evaluates_all_local_peaks(ctx->params->mode)))) {
                append_peak(scratch, ctx->params, peak_freq, peak_magnitude, ctx->df);
            }
        }
    }

    work->nPeaks = scratch->nPeaks;
    work->status = 0;
}

static inline bool execute_direct_grid(buffer_t *buffer, parameters *params, const eval_method_t *eval_method, kt_forpool_t *direct_pool, double fmin,
                                       double fstep, uint32_t nfreq, float threshold, bool threshold_valid, double df, int precision, bool scan_peaks,
                                       bool write_spectrum) {
    if (!buffer->power) return false;
    uint32_t chunk = direct_grid_chunk_size(buffer->n);
    uint32_t nworksets = chunk > 0U ? (nfreq + chunk - 1U) / chunk : 0U;
    if (nworksets == 0U) return true;
    int nthreads = (direct_pool && direct_pool->n_threads > 1) ? direct_pool->n_threads : 1;
    if (nthreads < 1) nthreads = 1;

    direct_grid_workset_t *worksets = calloc(nworksets, sizeof(*worksets));
    buffer_t *scratch = calloc((size_t)nthreads, sizeof(*scratch));
    if (!worksets || !scratch) {
        free(worksets);
        free(scratch);
        return false;
    }

    bool ok = true;
    for (int t = 0; t < nthreads; ++t) {
        if (alloc_direct_scratch(&scratch[t], buffer, params, false) != 0) {
            ok = false;
            break;
        }
    }

    int max_peaks = params->npeaks > 0 ? params->npeaks : 0;
    for (uint32_t w = 0; ok && w < nworksets; ++w) {
        uint32_t core_begin = w * chunk;
        uint32_t core_end = core_begin + chunk;
        if (core_end > nfreq) core_end = nfreq;
        uint32_t eval_begin = core_begin > 0U ? core_begin - 1U : core_begin;
        uint32_t eval_end = core_end < nfreq ? core_end + 1U : core_end;
        worksets[w].core_begin = core_begin;
        worksets[w].core_end = core_end;
        worksets[w].eval_begin = eval_begin;
        worksets[w].eval_end = eval_end;
        worksets[w].status = -1;
        worksets[w].spectrum = write_spectrum ? sdsempty() : NULL;
        worksets[w].values = malloc((size_t)(eval_end - eval_begin) * sizeof(float));
        worksets[w].peaks = max_peaks > 0 ? calloc((size_t)max_peaks, sizeof(peak_t)) : NULL;
        if (!worksets[w].values || (write_spectrum && !worksets[w].spectrum) || (max_peaks > 0 && !worksets[w].peaks)) ok = false;
    }

    if (ok) {
        direct_grid_context_t ctx = {.buffer = buffer,
                                     .params = params,
                                     .eval_method = eval_method,
                                     .threshold_valid = threshold_valid,
                                     .scan_peaks = scan_peaks,
                                     .write_spectrum = write_spectrum,
                                     .threshold = threshold,
                                     .fmin = fmin,
                                     .fstep = fstep,
                                     .df = df,
                                     .nfreq = nfreq,
                                     .precision = precision,
                                     .scratch = scratch,
                                     .worksets = worksets};
        if (direct_pool && nthreads > 1) {
            kt_forpool(direct_pool, direct_grid_worker, &ctx, (long)nworksets);
        } else {
            for (uint32_t w = 0; w < nworksets; ++w) direct_grid_worker(&ctx, (long)w, 0);
        }
        for (uint32_t w = 0; w < nworksets; ++w) {
            if (worksets[w].status != 0) ok = false;
        }
    }

    if (ok) {
        buffer->nPeaks = 0U;
        if (scan_peaks && max_peaks > 0) {
            for (uint32_t w = 0; w < nworksets; ++w) {
                for (uint32_t p = 0; p < worksets[w].nPeaks; ++p) {
                    direct_merge_peak_stable(buffer->peaks, &buffer->nPeaks, max_peaks, &worksets[w].peaks[p], params, eval_method);
                }
            }
        }
        if (write_spectrum) {
            sdsclear(buffer->spectrum);
            for (uint32_t w = 0; w < nworksets; ++w) {
                buffer->spectrum = sdscatsds(buffer->spectrum, worksets[w].spectrum);
            }
        }
    }

    for (uint32_t w = 0; w < nworksets; ++w) {
        sdsfree(worksets[w].spectrum);
        free(worksets[w].values);
        free(worksets[w].peaks);
    }
    for (int t = 0; t < nthreads; ++t) free_direct_scratch(&scratch[t]);
    free(worksets);
    free(scratch);
    return ok;
}

void process_target(char *in_file, buffer_t *buffer, parameters *params, const bool batch, kt_forpool_t *direct_pool) {
    PROFILE_COUNT_TARGET();
    PROFILE_START(read);
    if (has_fits_suffix(in_file)) {
        read_fits(in_file, buffer);
    } else {
        read_dat(in_file, buffer);
    }
    PROFILE_END(read, PROFILE_PHASE_READ);
    PROFILE_START(preprocess);
    preprocess_buffer(buffer, params->epsilon, params->mode);
    PROFILE_END(preprocess, PROFILE_PHASE_PREPROCESS);

    peak_t *peak_base = buffer->peaks;
    peak_t *selected_peaks = peak_base;
    peak_t *candidate_peaks = params->prewhiten ? peak_base + params->npeaks : peak_base;
    int selected_count = 0;

    const int n = 1 + (int)(log10(buffer->x[buffer->n - 1] * (double)(params->oversamplingFactor * params->nterms)));
    double fmin = effective_grid_fmin(params, buffer->x[buffer->n - 1]);
    const bool use_aov = periodogram_uses_aov(params->periodogramMethod);
    const eval_method_t *eval_method = eval_method_for_params(params);
    const bool direct_eval_grid = mode_uses_direct_eval_grid(params->mode);
    const int aov_n_eff = periodogram_effective_n(buffer);
    const bool aov_valid = !use_aov || aov_target_has_dof(buffer, params->nterms);
    const bool threshold_valid = direct_eval_grid || aov_valid;
    const float threshold = threshold_valid ? (float)correct_threshold(params, buffer) : INFINITY;
    double fstep = 0.0;
    if (mode_uses_direct_eval_grid(params->mode)) {
        fstep = eval_method->direct_frequency_step(buffer->n, buffer->x[buffer->n - 1], params);
    } else {
        fstep = 1.0 / (double)(params->nterms * (double)params->oversamplingFactor * buffer->x[buffer->n - 1] * 0.5);
    }
    double df = 1.0 / buffer->x[buffer->n - 1];
    uint32_t nfreq = target_frequency_count(params, fmin, fstep);
    if (nfreq > buffer->maxFreqCount) nfreq = buffer->maxFreqCount;
    if (!mode_uses_direct_gb_grid(params->mode) && !activate_target_nufft_plan(buffer, params, nfreq)) {
        fprintf(stderr, "Failed to activate NuFFT plan for %s\n", in_file);
        goto cleanup;
    }
    const bool spectrum_prewhiten = params->spectrum && params->prewhiten;
    float *spectrum_matrix = NULL;
    uint32_t spectrum_columns = 0U;
    uint32_t spectrum_capacity = 0U;
    if (spectrum_prewhiten) {
        int requested_peaks = params->npeaks > 0 ? params->npeaks : 0;
        spectrum_capacity = (uint32_t)requested_peaks + 1U;
        if (nfreq == 0U) {
            spectrum_matrix = malloc(sizeof(float));
        } else if (spectrum_capacity > 0U && (size_t)nfreq <= SIZE_MAX / (size_t)spectrum_capacity / sizeof(float)) {
            spectrum_matrix = malloc((size_t)nfreq * (size_t)spectrum_capacity * sizeof(float));
        }
        if (!spectrum_matrix) {
            fprintf(stderr, "Failed to allocate prewhitened spectrum matrix for %s\n", in_file);
            goto cleanup;
        }
    }

    if (!batch && params->debug) {
        printf("\tNumber of target frequencies: %u\n\n", nfreq);
    }

    char stringBuff[64];
    bool aov_streamed = false;

    do {
        aov_streamed = false;
        buffer->peaks = candidate_peaks;
        buffer->nPeaks = 0;
        bool direct_grid_evaluated = false;
        if (mode_uses_direct_eval_grid(params->mode)) {
            direct_grid_evaluated = true;
            PROFILE_START(peak_scan);
            if (!execute_direct_grid(buffer, params, eval_method, direct_pool, fmin, fstep, nfreq, threshold, threshold_valid, df, n, true,
                                     params->spectrum && !spectrum_prewhiten)) {
                fprintf(stderr, "Direct evaluation grid failed for %s\n", in_file);
            }
            PROFILE_END(peak_scan, PROFILE_PHASE_PEAK_SCAN);
        } else if (!mode_uses_direct_gb_grid(params->mode)) {
            if (use_aov) {
                if (direct_pool && direct_pool->n_threads > 1 && params->buffers && params->nbuffers >= direct_pool->n_threads) {
                    // Parallel AoV: fills power[], serial post-pass handles spectrum + peaks
                    if (!execute_aov_sweep_parallel(buffer, params, direct_pool, fmin, fstep, nfreq)) {
                        fprintf(stderr, "AoV periodogram failed for %s\n", in_file);
                    }
                } else {
                    aov_streamed = true;
                    if (execute_aov_sweep(buffer, params, fmin, fstep, nfreq, threshold, df, params->spectrum && !spectrum_prewhiten, params->spectrum, true, n,
                                          stringBuff) != 0) {
                        fprintf(stderr, "AoV periodogram failed for %s\n", in_file);
                    }
                }
            } else {
                if (direct_pool && direct_pool->n_threads > 1 && params->buffers && params->nbuffers >= direct_pool->n_threads) {
                    // Parallel NuFFT: fills power[], serial post-pass handles spectrum + peaks
                    if (!execute_nufft_sweep_parallel(buffer, params, direct_pool, fmin, fstep, nfreq)) {
                        fprintf(stderr, "IHS periodogram failed for %s\n", in_file);
                    }
                } else {
                    execute_nufft_sweep(buffer, params, fmin, fstep, nfreq);
                }
            }
        }

        if (spectrum_prewhiten) {
            PROFILE_START(spectrum);
            if (spectrum_columns < spectrum_capacity) {
                capture_spectrum_column(buffer, params, eval_method, use_aov, direct_eval_grid, aov_n_eff, fmin, fstep, nfreq, params->nterms, spectrum_matrix,
                                        spectrum_columns++);
            }
            PROFILE_END(spectrum, PROFILE_PHASE_PEAK_SCAN);
        } else if (params->spectrum && !aov_streamed && !direct_grid_evaluated) {
            PROFILE_START(spectrum);
            for (uint32_t i = 0; i < nfreq; ++i) {
                double freq = fmin + ((double)i * fstep);
                float magnitude;
                if (direct_eval_grid) {
                    magnitude = get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL);
                } else if (use_aov) {
                    float r2 = mode_uses_direct_eval_grid(params->mode) ? aov_get_stat(buffer, params, freq) : buffer->power[i];
                    magnitude = aov_likelihood_from_r2(r2, params->nterms, aov_n_eff);
                } else {
                    magnitude = mode_uses_direct_eval_grid(params->mode) ? get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL)
                                                                         : correct_ihs_res(buffer->power[i], params->nterms);
                }
                if (!float_is_finite_bits(magnitude)) magnitude = 0.0f;
                appendFreq(freq, magnitude, n, params->outputPeriod, &buffer->spectrum, &stringBuff[0]);
            }
            PROFILE_END(spectrum, PROFILE_PHASE_PEAK_SCAN);
        }

        if (!aov_streamed && !direct_grid_evaluated) {
            PROFILE_START(peak_scan);
            for (uint32_t i = 1; i + 1 < nfreq; ++i) {
                double freq = fmin + ((double)i * fstep);
                float magnitude = 0.0f;
                float left = 0.0f;
                float right = 0.0f;
                if (!mode_uses_direct_gb_grid(params->mode)) {
                    magnitude = buffer->power[i];
                    left = buffer->power[i - 1];
                    right = buffer->power[i + 1];
                }
                if (direct_eval_grid) {
                    magnitude = get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL);
                    left = get_eval_likelihood(buffer, params, eval_method, freq - fstep, NULL, NULL);
                    right = get_eval_likelihood(buffer, params, eval_method, freq + fstep, NULL, NULL);
                } else if (use_aov && mode_uses_direct_eval_grid(params->mode)) {
                    magnitude = aov_get_stat(buffer, params, freq);
                    left = aov_get_stat(buffer, params, freq - fstep);
                    right = aov_get_stat(buffer, params, freq + fstep);
                } else if (mode_uses_direct_eval_grid(params->mode)) {
                    magnitude = get_eval_likelihood(buffer, params, eval_method, freq, NULL, NULL);
                    left = get_eval_likelihood(buffer, params, eval_method, freq - fstep, NULL, NULL);
                    right = get_eval_likelihood(buffer, params, eval_method, freq + fstep, NULL, NULL);
                }
                double peak_freq = freq;
                float peak_magnitude = magnitude;
                if (!quadratic_peak_position(freq, fstep, left, magnitude, right, params->oversamplingFactor, &peak_freq, &peak_magnitude)) continue;
                if (magnitude > left && magnitude > right && (magnitude > threshold || (threshold_valid && mode_evaluates_all_local_peaks(params->mode)))) {
                    if (use_aov && !direct_eval_grid) {
                        aov_append_peak(buffer, params, peak_freq, peak_magnitude, df);
                    } else {
                        append_peak(buffer, params, peak_freq, peak_magnitude, df);
                    }
                }
            }
            PROFILE_END(peak_scan, PROFILE_PHASE_PEAK_SCAN);
        }

        int candidate_count = (int)buffer->nPeaks;
        PROFILE_START(sort);
        if (use_aov && !direct_eval_grid) {
            aov_sort_peaks(buffer->peaks, candidate_count, buffer, params, df);
        } else {
            sortPeaks(buffer->peaks, candidate_count, buffer, params, df);
        }
        PROFILE_END(sort, PROFILE_PHASE_SORT);

        if (!params->prewhiten) break;
        if (candidate_count <= 0 || selected_count >= params->npeaks) break;

        peak_t selected = buffer->peaks[0];
        if (selected.p <= 0.0f || !float_is_finite_bits(selected.p)) break;

        eval_result_t prewhiten_result = eval_prewhiten(buffer, params, eval_method, selected.freq);
        float threshold_stat = prewhiten_result.stat;
        if (!float_is_finite_bits(threshold_stat) || threshold_stat < params->r2_threshold) break;

        selected_peaks[selected_count++] = selected;
        refresh_weighted_signal_buffer(buffer, params->epsilon);
    } while (params->prewhiten && selected_count < params->npeaks);

    if (spectrum_prewhiten && selected_count > 0 && spectrum_columns == (uint32_t)selected_count && spectrum_columns < spectrum_capacity) {
        bool final_aov_streamed = false;
        if (mode_uses_direct_eval_grid(params->mode)) {
            if (!execute_direct_grid(buffer, params, eval_method, direct_pool, fmin, fstep, nfreq, threshold, threshold_valid, df, n, false, false)) {
                fprintf(stderr, "Direct evaluation grid failed for %s\n", in_file);
            }
        } else if (!mode_uses_direct_gb_grid(params->mode)) {
            if (use_aov) {
                if (direct_pool && direct_pool->n_threads > 1 && params->buffers && params->nbuffers >= direct_pool->n_threads) {
                    if (!execute_aov_sweep_parallel(buffer, params, direct_pool, fmin, fstep, nfreq)) {
                        fprintf(stderr, "AoV periodogram failed for %s\n", in_file);
                    }
                } else {
                    final_aov_streamed = true;
                    if (execute_aov_sweep(buffer, params, fmin, fstep, nfreq, threshold, df, false, true, false, n, stringBuff) != 0) {
                        fprintf(stderr, "AoV periodogram failed for %s\n", in_file);
                    }
                }
            } else {
                if (direct_pool && direct_pool->n_threads > 1 && params->buffers && params->nbuffers >= direct_pool->n_threads) {
                    if (!execute_nufft_sweep_parallel(buffer, params, direct_pool, fmin, fstep, nfreq)) {
                        fprintf(stderr, "IHS periodogram failed for %s\n", in_file);
                    }
                } else {
                    execute_nufft_sweep(buffer, params, fmin, fstep, nfreq);
                }
            }
        }
        (void)final_aov_streamed;
        capture_spectrum_column(buffer, params, eval_method, use_aov, direct_eval_grid, aov_n_eff, fmin, fstep, nfreq, params->nterms, spectrum_matrix,
                                spectrum_columns++);
    }

    if (params->prewhiten) {
        buffer->peaks = selected_peaks;
        buffer->nPeaks = (uint32_t)selected_count;
    }

    {
        PROFILE_START(output);
        if (!batch) {
            print_peaks(buffer, params, n, stringBuff, in_file, params->mode, params->gbEvalMode);
        } else {
            append_peaks(buffer, params, n, stringBuff, in_file, params->mode, params->gbEvalMode);
        }
        if (params->spectrum) {
            if (spectrum_prewhiten) {
                write_spectrum_matrix_tsv(in_file, fmin, fstep, nfreq, spectrum_columns, spectrum_matrix, n, params->outputPeriod, stringBuff);
            } else {
                write_tsv(buffer, in_file);
            }
        }
        PROFILE_END(output, PROFILE_PHASE_OUTPUT);
    }
    free(spectrum_matrix);
    buffer->peaks = peak_base;
    buffer->nPeaks = 0;
    PROFILE_FLUSH_THREAD();
    return;

cleanup:
    free(spectrum_matrix);
    buffer->peaks = peak_base;
    buffer->nPeaks = 0;
    PROFILE_FLUSH_THREAD();
}

void process_targets(void *data, long i, int thread_id) {
    parameters *params = (parameters *)data;

    if (!params->buffers[thread_id]->allocated) {
        alloc_buffer(params->buffers[thread_id], params);
    }
    process_target(kv_A(params->targets, i).path, params->buffers[thread_id], params, true, NULL);

    int total = (int)kv_size(params->targets);
    int permile = total / 1000;
    if (permile == 0) permile = 1;
    int current = __sync_add_and_fetch(&params->iter_count, 1);
    if (current % permile == 0 || current >= total) {
        float progress = (float)current * 100.0f / (float)total;
        double elapsed = elapsed_seconds_since(&params->batch_start_time);
        char time_buf[64];
        if (elapsed > 0.001 && current > 0) {
            double remaining_sec = (double)(total - current) * elapsed / (double)current;
            format_time_remaining(remaining_sec, time_buf, sizeof(time_buf));
        } else {
            snprintf(time_buf, sizeof(time_buf), "calculating...");
        }
        printf("\rComputation in progress: %.1f%% complete | %s%s", progress, time_buf, current >= total ? "\n" : "");
        fflush(stdout);
    }
}

#endif
