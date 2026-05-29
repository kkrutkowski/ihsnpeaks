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

static inline void fill_twiddle_ladder(buffer_t *buffer, const parameters *params, double fstep, double freq_factor) {
    for (uint32_t level = 0; level < params->ladderLevels; ++level) {
        size_t offset = (size_t)level * (size_t)buffer->paddedLen;
        double advance = nufft1_twiddle_ladder_advance((int)params->outputLen, (int)level);
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

static inline void reset_work_ladder(buffer_t *buffer, const parameters *params) {
    for (uint32_t level = 0; level < params->ladderLevels; ++level) {
        size_t offset = (size_t)level * (size_t)buffer->paddedLen;
        memcpy(buffer->workReal + offset, buffer->inputReal, (size_t)buffer->paddedLen * sizeof(float));
        memcpy(buffer->workImag + offset, buffer->inputImag, (size_t)buffer->paddedLen * sizeof(float));
    }
}

static inline void advance_work_ladder(buffer_t *buffer, const parameters *params, size_t next_block) {
    int level = nufft1_twiddle_ladder_carry_level(next_block, (int)params->ladderLevels);
    if (level >= (int)params->ladderLevels) level = (int)params->ladderLevels - 1;
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
    memset(buffer->power, 0, ((size_t)nfreq + 1U) * sizeof(float));
    nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);

    size_t num_blocks = ((size_t)nfreq + (size_t)params->outputLen - 1U) / (size_t)params->outputLen;
    for (int t = 0; t < buffer->terms; ++t) {
        double freq_factor = (double)(t + 1);
        fill_complex_input(buffer, fmin, freq_factor);
        fill_twiddle_ladder(buffer, params, fstep, freq_factor);
        reset_work_ladder(buffer, params);

        for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
            uint32_t base = (uint32_t)(block_idx * (size_t)params->outputLen);
            uint32_t count = params->outputLen;
            if (base + count > nfreq) count = nfreq - base;
            nufft1_execute(buffer->nufftWorkspace, buffer->workReal, buffer->workImag, buffer->blockReal, buffer->blockImag, t + 1);
            accumulate_power(buffer, base, count);
            if (block_idx + 1U < num_blocks) {
                advance_work_ladder(buffer, params, block_idx + 1U);
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
    read_dat(in_file, buffer);
    preprocess_buffer(buffer, params->epsilon, params->mode);

    int prewhitening_iter = 0;
    double *backup_r2 = calloc(params->npeaks, sizeof(double));
    double *backup_p = calloc(params->npeaks, sizeof(double));

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

    if (!batch && params->debug) {
        printf("\tNumber of target frequencies: %u\n\n", nfreq);
    }

    char stringBuff[32];
    bool aov_streamed = false;

prewhiten:
    aov_streamed = false;
    if (prewhitening_iter > 0) {
        double prewhitening_freq = buffer->peaks[prewhitening_iter - 1].freq;
        float prewhiten_r2 = get_r2(buffer, prewhitening_freq, NULL, true, params->gbAlpha);

        backup_r2[prewhitening_iter] = buffer->peaks[prewhitening_iter - 1].r2;
        backup_p[prewhitening_iter] = buffer->peaks[prewhitening_iter - 1].p;

        double ysum = 0;
        for (unsigned int i = 0; i < buffer->n; i++) {
            ysum += fabs(buffer->y[i]);
        }
        double norm_neff = buffer->amp_neff > 0.0 ? buffer->amp_neff : buffer->neff;
        float norm = (sqrtf((float)norm_neff) / (float)(ysum));
        for (unsigned int i = 0; i < buffer->n; i++) {
            buffer->wy[i] = buffer->y[i] * norm;
        }

        float threshold_stat = params->gbEvalMode == GB_EVAL_GBAS ? prewhiten_r2 : buffer->peaks[prewhitening_iter - 1].r2;
        if (threshold_stat < params->r2_threshold) {
            buffer->nPeaks -= 1;
            goto next;
        }
    }

    buffer->nPeaks = prewhitening_iter;
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
            appendFreq(freq, magnitude, n, &buffer->spectrum, &stringBuff[0]);
        }
    }

    if (!aov_streamed)
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

    if (use_aov) {
        aov_sort_peaks(&buffer->peaks[prewhitening_iter], buffer->nPeaks, buffer, params, df);
    } else {
        sortPeaks(&buffer->peaks[prewhitening_iter], buffer->nPeaks, buffer, params->mode, df, params->nterms, params->gbEvalMode, params->gbAlpha);
    }

    if (params->prewhiten && prewhitening_iter < buffer->nPeaks) {
        prewhitening_iter += 1;
        buffer->nPeaks = prewhitening_iter;
        goto prewhiten;
    }

next:
    if (prewhitening_iter > 1) {
        for (int i = 1; i <= prewhitening_iter; i++) {
            buffer->peaks[i - 1].r2 = backup_r2[i];
            buffer->peaks[i - 1].p = backup_p[i];
        }
        buffer->nPeaks = prewhitening_iter;
    }

    if (!batch) {
        print_peaks(buffer, params, n, stringBuff, in_file, params->mode, params->gbEvalMode);
    } else {
        append_peaks(buffer, params, n, stringBuff, in_file, params->mode, params->gbEvalMode);
    }
    if (params->spectrum) {
        write_tsv(buffer, in_file);
    }
    free(backup_r2);
    free(backup_p);
    buffer->nPeaks = 0;
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
