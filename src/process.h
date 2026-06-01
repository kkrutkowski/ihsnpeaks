#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "profile.h"
#include "utils/readout.h"
#include "utils/simd.h"
#include "utils/strings.h"

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

#include "utils/aov.h"

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

void process_target(char *in_file, buffer_t *buffer, parameters *params, const bool batch) {
    PROFILE_COUNT_TARGET();
    PROFILE_START(read);
    read_dat(in_file, buffer);
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
    const int aov_n_eff = periodogram_effective_n(buffer);
    const bool aov_valid = !use_aov || aov_target_has_dof(buffer, params->nterms);
    const float threshold = aov_valid ? (float)correct_threshold(params, buffer) : INFINITY;
    double fstep = mode_uses_direct_gb_grid(params->mode)
                       ? gb_direct_frequency_step(buffer->n, buffer->x[buffer->n - 1], params->oversamplingFactor, params->gbAlpha)
                       : 1.0 / (double)(params->nterms * (double)params->oversamplingFactor * buffer->x[buffer->n - 1] * 0.5);
    double df = 1.0 / buffer->x[buffer->n - 1];
    uint32_t nfreq = target_frequency_count(params, fmin, fstep);
    if (nfreq > buffer->maxFreqCount) nfreq = buffer->maxFreqCount;
    if (!mode_uses_direct_gb_grid(params->mode) && !activate_target_nufft_plan(buffer, params, nfreq)) {
        fprintf(stderr, "Failed to activate NuFFT plan for %s\n", in_file);
        goto cleanup;
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
        if (!mode_uses_direct_gb_grid(params->mode)) {
            if (use_aov) {
                aov_streamed = true;
                if (execute_aov_sweep(buffer, params, fmin, fstep, nfreq, threshold, df, params->spectrum, n, stringBuff) != 0) {
                    fprintf(stderr, "AoV periodogram failed for %s\n", in_file);
                }
            } else {
                execute_nufft_sweep(buffer, params, fmin, fstep, nfreq);
            }
        }

        if (params->spectrum && !aov_streamed) {
            PROFILE_START(spectrum);
            for (uint32_t i = 0; i < nfreq; ++i) {
                double freq = fmin + ((double)i * fstep);
                float magnitude;
                if (use_aov) {
                    float r2 = mode_uses_direct_gb_grid(params->mode) ? aov_get_stat(buffer, params, freq) : buffer->power[i];
                    magnitude = aov_likelihood_from_r2(r2, params->nterms, aov_n_eff);
                } else {
                    magnitude = mode_uses_direct_gb_grid(params->mode) ? get_gb_likelihood(buffer, freq, NULL, NULL, params->gbEvalMode, params->gbAlpha)
                                                                       : correct_ihs_res(buffer->power[i], params->nterms);
                }
                if (!float_is_finite_bits(magnitude)) magnitude = 0.0f;
                appendFreq(freq, magnitude, n, params->outputPeriod, &buffer->spectrum, &stringBuff[0]);
            }
            PROFILE_END(spectrum, PROFILE_PHASE_PEAK_SCAN);
        }

        if (!aov_streamed) {
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
                if (use_aov && mode_uses_direct_gb_grid(params->mode)) {
                    magnitude = aov_get_stat(buffer, params, freq);
                    left = aov_get_stat(buffer, params, freq - fstep);
                    right = aov_get_stat(buffer, params, freq + fstep);
                } else if (mode_uses_direct_gb_grid(params->mode)) {
                    magnitude = get_gb_likelihood(buffer, freq, NULL, NULL, params->gbEvalMode, params->gbAlpha);
                    left = get_gb_likelihood(buffer, freq - fstep, NULL, NULL, params->gbEvalMode, params->gbAlpha);
                    right = get_gb_likelihood(buffer, freq + fstep, NULL, NULL, params->gbEvalMode, params->gbAlpha);
                }
                double peak_freq = freq;
                float peak_magnitude = magnitude;
                if (!quadratic_peak_position(freq, fstep, left, magnitude, right, params->oversamplingFactor, &peak_freq, &peak_magnitude)) continue;
                if (magnitude > left && magnitude > right && (magnitude > threshold || (aov_valid && mode_evaluates_all_local_peaks(params->mode)))) {
                    if (use_aov) {
                        aov_append_peak(buffer, params, peak_freq, peak_magnitude, df);
                    } else {
                        append_peak(buffer, params->npeaks, params->mode, peak_freq, peak_magnitude, df, params->gbEvalMode, params->gbAlpha);
                    }
                }
            }
            PROFILE_END(peak_scan, PROFILE_PHASE_PEAK_SCAN);
        }

        int candidate_count = (int)buffer->nPeaks;
        PROFILE_START(sort);
        if (use_aov) {
            aov_sort_peaks(buffer->peaks, candidate_count, buffer, params, df);
        } else {
            sortPeaks(buffer->peaks, candidate_count, buffer, params->mode, df, params->nterms, params->gbEvalMode, params->gbAlpha);
        }
        PROFILE_END(sort, PROFILE_PHASE_SORT);

        if (!params->prewhiten) break;
        if (candidate_count <= 0 || selected_count >= params->npeaks) break;

        peak_t selected = buffer->peaks[0];
        if (selected.p <= 0.0f || !float_is_finite_bits(selected.p)) break;

        float prewhiten_r2 = get_r2(buffer, selected.freq, NULL, true, params->gbAlpha);
        float threshold_stat = (params->gbEvalMode == GB_EVAL_GBAW || params->mode == 0) ? prewhiten_r2 : selected.r2;
        if (!float_is_finite_bits(threshold_stat) || threshold_stat < params->r2_threshold) break;

        selected_peaks[selected_count++] = selected;
        refresh_weighted_signal_buffer(buffer, params->epsilon);
    } while (params->prewhiten && selected_count < params->npeaks);

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
            write_tsv(buffer, in_file);
        }
        PROFILE_END(output, PROFILE_PHASE_OUTPUT);
    }
    buffer->peaks = peak_base;
    buffer->nPeaks = 0;
    PROFILE_FLUSH_THREAD();
    return;

cleanup:
    buffer->peaks = peak_base;
    buffer->nPeaks = 0;
    PROFILE_FLUSH_THREAD();
}

void process_targets(void *data, long i, int thread_id) {
    parameters *params = (parameters *)data;

    int permile = kv_size(params->targets) / 1000;
    if (permile == 0) {
        permile += 1;
    }
    pthread_mutex_lock(&params->counter_mutex);
    params->iter_count += 1;
    int current_iter = params->iter_count;
    if (current_iter % permile == 0 || current_iter == kv_size(params->targets)) {
        float progress = (float)(current_iter) * 100.0f / (float)kv_size(params->targets);
        printf("Computation in progress: %.1f%% complete\r", progress);
        fflush(stdout);
    }
    pthread_mutex_unlock(&params->counter_mutex);
    if (!params->buffers[thread_id]->allocated) {
        alloc_buffer(params->buffers[thread_id], params);
    }
    process_target(kv_A(params->targets, i).path, params->buffers[thread_id], params, true);
}

#endif
