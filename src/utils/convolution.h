#ifndef FTEST_H

#    include <fdist.h>
#    include <float.h>
#    include <math.h>
#    include <stddef.h>
#    include <stdint.h>
#    include <stdio.h>
#    include <stdlib.h>
#    include <string.h>

#    include "common.h"
#    include "simd.h"

// Computes the values of self-correlation of the convolved array elements
static inline double corr(int r) {
    double L = (2.0 * (double)r) + 1.0;
    return (((2 * L * L) + 1) / (3 * L * L * L));
}

static inline double clamp_gbls_scale(double scale) {
    if (scale <= 0.0) {
        return 0.0;
    }
    double range = 4.0;

    // Brute-force approach. Not used since sharp edges negatively affect final interpolation step.
    // if (scale < 1.0 / range) return 1.0 / range; if (scale > range) return range;
    // return scale;

    // let's use tanh() here instead
    double logrange = log2(range);
    double arg = log2(scale) / logrange;
    return pow(2.0, tanh(arg) * logrange);
}

static inline int32_t wrapidx(int32_t idx, int32_t n) {
    while (idx < 0) {
        idx += n;
    }
    while (idx >= n) {
        idx -= n;
    }
    return idx;
}

// Extract key from the 64-bit value (16-bit key at the lower 16 bits)
static inline uint16_t extract_key(uint64_t value) { return (uint16_t)value; }  // Assuming lowest 16 bits only on little endian

static inline uint16_t phase_key_10(double phase) { return (uint16_t)(((uint16_t)(phase * 1024.0)) & 1023U); }

static inline void csort64_10(uint64_t** array, size_t n, uint64_t** aux_buffer, size_t* indices) {
    // Clear the indices array (initialized to 0)
    memset(indices, 0, 1024 * sizeof(size_t));

    // Compute the bin sizes (counting phase)
    for (size_t i = 0; i < n; i++) {
        size_t index_tmp = extract_key((*array)[i]);
        indices[index_tmp] += 1;
    }

    // Compute the first index of each bin (prefix sum phase)
    for (size_t i = 1; i < 1024; i++) {
        indices[i] += indices[i - 1];
    }

    // Rearrange elements into the aux_buffer based on indices (placement phase)
    for (size_t i = 0; i < n; i++) {
        uint16_t index_tmp = extract_key((*array)[i]);
        size_t pos = --indices[index_tmp];  // Decrease the position before placing the element
        (*aux_buffer)[pos] = (*array)[i];
    }
}

void convolve(kvpair* in, double* temp, double* out, int r, int n) {
    for (int j = 0; j < n; j++) {
        out[j] = in[j].parts.val;
    }
    double width = (2.0 * (double)r) + 1.0;
    double inv_width = 1.0 / width;
    double norm = inv_width * inv_width;
    norm *= norm;

    for (int i = 0; i <= 3; i++) {
        double val = 0.0;
        int idx_hi = r;
        int idx_lo = n - r;
        if (i % 2 == 0) {
            for (int j = n - (r); j < n + (r); j++) {
                val += out[wrapidx(j, n)];
            }
            for (int j = 0; j < n; j++) {
                if (idx_hi >= n) {
                    idx_hi = wrapidx(idx_hi, n);
                }
                if (idx_lo >= n) {
                    idx_lo = wrapidx(idx_lo, n);
                }
                val += out[idx_hi];
                temp[j] = val;
                val -= out[idx_lo];
                idx_hi += 1;
                idx_lo += 1;
            }
        } else {
            for (int j = n - (r); j < n + (r); j++) {
                val += temp[wrapidx(j, n)];
            }
            for (int j = 0; j < n; j++) {
                if (idx_hi >= n) {
                    idx_hi = wrapidx(idx_hi, n);
                }
                if (idx_lo >= n) {
                    idx_lo = wrapidx(idx_lo, n);
                }
                val += temp[idx_hi];
                out[j] = val;
                val -= temp[idx_lo];
                idx_hi += 1;
                idx_lo += 1;
            }
        }
    }
    for (int j = 0; j < n; j++) {
        out[j] *= norm;
    }  // Normalize
}

typedef struct {
    kvpair* sorted;
    double* model;
    double scale;
    double min_model;
    double max_model;
    double Sxx;
    double Syy;
    double Sxy;
    int n;
} gb_projection_t;

static inline bool gb_prepare_projection(buffer_t* buffer, double freq, bool keep_indices, float gbAlpha, double model_multiplier,
                                         gb_projection_t* projection) {
    if (!buffer || !projection || !buffer->buf || !buffer->buf[0] || !buffer->buf[1] || !buffer->buf[2] || !buffer->pidx || buffer->n == 0) return false;

    double min_model = 0.0;
    double max_model = 0.0;
    double Sxx = 0.0;
    double Syy = 0.0;
    double Sxy = 0.0;

    kvpair* input = (kvpair*)(buffer->buf[0]);
    kvpair* sorted = (kvpair*)(buffer->buf[1]);
    double* tmp = (double*)(buffer->buf[0]);
    double* output = (double*)(buffer->buf[2]);
    uint64_t** sort_in = (uint64_t**)(&buffer->buf[0]);
    uint64_t** sort_out = (uint64_t**)(&buffer->buf[1]);

    for (int i = 0; i < buffer->n; i++) {
        input[i].parts.key = phase_key_10(buffer->x[i] * freq);
        input[i].parts.val = buffer->y[i];
        if (keep_indices) {
            input[i].parts.idx = i;
        }
    }

    csort64_10(sort_in, buffer->n, sort_out, buffer->pidx);  // sort the pairs by phase
    int r = gb_convolution_radius(buffer->n, gbAlpha);

    convolve(sorted, tmp, output, r, buffer->n);

    double correction = corr(r);
    for (int i = 0; i < buffer->n; i++) {
        double y_val = (double)(sorted[i].parts.val);
        double smooth_val = (output[i] - (correction * y_val)) * model_multiplier;
        output[i] = smooth_val;
        if (i == 0 || smooth_val < min_model) min_model = smooth_val;
        if (i == 0 || smooth_val > max_model) max_model = smooth_val;

        Sxx += smooth_val * smooth_val;
        Syy += y_val * y_val;
        Sxy += y_val * smooth_val;
    }

    double scale = 0.0;
    if (Sxx > 0.0 && Syy > 0.0) {
        scale = clamp_gbls_scale(Sxy / Sxx);
    }

    projection->sorted = sorted;
    projection->model = output;
    projection->scale = scale;
    projection->min_model = min_model;
    projection->max_model = max_model;
    projection->Sxx = Sxx;
    projection->Syy = Syy;
    projection->Sxy = Sxy;
    projection->n = (int)buffer->n;
    return true;
}

static inline float get_r2(buffer_t* buffer, double freq, float* amp, bool prewhiten, float gbAlpha) {
    gb_projection_t projection = {0};
    if (!gb_prepare_projection(buffer, freq, prewhiten, gbAlpha, 1.0, &projection)) return 0.0f;

    float R2 = 0.0f;
    if (projection.Sxx > 0.0 && projection.Syy > 0.0) {
        double scale = projection.scale;
        double Sxy = projection.Sxy;
        double Sxx = projection.Sxx;
        double Syy = projection.Syy;
        double explained = ((2.0 * scale * Sxy) - (scale * scale * Sxx)) / Syy;
        R2 = (float)explained;
        if (R2 < 0.0f) R2 = 0.0f;
        if (R2 > (1.0f - 1e-15f)) R2 = 1.0f - 1e-15f;
    }
    if (amp) {
        *amp = (float)(fabs(projection.scale) * (projection.max_model - projection.min_model));
    }
    if (prewhiten) {
        for (int i = 0; i < projection.n; i++) {
            buffer->y[projection.sorted[i].parts.idx] -= (float)(projection.scale * projection.model[i]);
        }
    }

    return R2;
}

static inline float get_gbaw_t(buffer_t* buffer, double freq, float* amp, float gbAlpha) {
    if (!buffer || buffer->n == 0) return 0.0f;
    double ref = 0;
    double res = 0;

    int r = gb_convolution_radius(buffer->n, gbAlpha);
    double correction = corr(r);
    double multiplier = 1.0 / (1.0 - correction);
    gb_projection_t projection = {0};
    if (!gb_prepare_projection(buffer, freq, false, gbAlpha, multiplier, &projection)) return 0.0f;

    for (int i = 0; i < projection.n; i++) {
        double y_val = (double)(projection.sorted[i].parts.val);
        double model = projection.model[i];
        res += (y_val - model) * (y_val - model);
        ref += (y_val + model) * (y_val + model);
    }

    if (amp) {
        *amp = (float)(fabs(projection.scale) * (projection.max_model - projection.min_model));
    }

    return (float)(ref / res);
}

typedef struct {
    float likelihood;
    float stat;
    float amp;
    float width;
    bool valid;
} eval_result_t;

typedef struct {
    double y;
    double y2;
} bls_sums_t;

static inline double bls_relative_width_at(const parameters* params, int idx) {
    double min_width = params->blsMinRelWidth > 0.0 ? params->blsMinRelWidth : 0.01;
    double max_width = params->blsMaxRelWidth > 0.0 ? params->blsMaxRelWidth : 0.5;
    int count = params->blsWidthCount > 0 ? params->blsWidthCount : 10;
    if (max_width < min_width) {
        double tmp = min_width;
        min_width = max_width;
        max_width = tmp;
    }
    if (count <= 1 || min_width == max_width) return min_width;
    double t = (double)idx / (double)(count - 1);
    return exp(log(min_width) + (log(max_width) - log(min_width)) * t);
}

static inline uint64_t bls_pack_phase_index(uint16_t key, uint32_t idx) { return (((uint64_t)idx) << 16U) | (uint64_t)(key & 1023U); }

static inline uint32_t bls_unpack_index(uint64_t value) { return (uint32_t)(value >> 16U); }

static inline void bls_add_index(const buffer_t* buffer, uint32_t idx, double sign, bls_sums_t* sums) {
    double y = (double)buffer->y[idx];
    sums->y += sign * y;
    sums->y2 += sign * y * y;
}

static inline void bls_add_sorted(const buffer_t* buffer, const uint64_t* sorted, int pos, double sign, bls_sums_t* sums) {
    int n = (int)buffer->n;
    if (pos >= n) pos %= n;
    if (pos < 0) pos = wrapidx(pos, n);
    bls_add_index(buffer, bls_unpack_index(sorted[pos]), sign, sums);
}

static inline float bls_float_value(double value) {
    if (!double_is_finite_bits(value) || value <= 0.0) return 0.0f;
    if (value > (double)FLT_MAX) return FLT_MAX;
    return (float)value;
}

static inline eval_result_t get_bls_result(buffer_t* buffer, const parameters* params, double freq, bool prewhiten) {
    eval_result_t result = {0};
    if (!buffer || !params || !buffer->buf || !buffer->buf[0] || !buffer->buf[1] || !buffer->pidx || buffer->n <= 2U) return result;

    int n = (int)buffer->n;
    uint64_t* input = (uint64_t*)buffer->buf[0];
    uint64_t* sorted = (uint64_t*)buffer->buf[1];
    uint64_t** sort_in = (uint64_t**)(&buffer->buf[0]);
    uint64_t** sort_out = (uint64_t**)(&buffer->buf[1]);
    bls_sums_t total = {0};
    for (int i = 0; i < n; ++i) {
        uint16_t key = phase_key_10(buffer->x[i] * freq);
        input[i] = bls_pack_phase_index(key, (uint32_t)i);
        bls_add_index(buffer, (uint32_t)i, 1.0, &total);
    }
    csort64_10(sort_in, buffer->n, sort_out, buffer->pidx);

    double ref_sse = total.y2 - ((total.y * total.y) / (double)n);
    if (!(ref_sse > 0.0) || !double_is_finite_bits(ref_sse)) return result;

    int width_count = params->blsWidthCount > 0 ? params->blsWidthCount : 10;
    double best_score = -1.0;
    double best_r2 = 0.0;
    double best_amp = 0.0;
    double best_in_level = 0.0;
    double best_out_level = 0.0;
    int best_start = -1;
    int best_width = 0;

    for (int width_idx = 0; width_idx < width_count; ++width_idx) {
        double rel_width = bls_relative_width_at(params, width_idx);
        if (!(rel_width > 0.0) || !double_is_finite_bits(rel_width)) continue;
        int width = (int)ceil(rel_width * (double)n);
        if (width < 1) width = 1;
        if (width >= n) width = n - 1;

        bls_sums_t in = {0};
        for (int j = 0; j < width; ++j) {
            bls_add_sorted(buffer, sorted, j, 1.0, &in);
        }

        int out_count = n - width;
        for (int start = 0; start < n; ++start) {
            bls_sums_t out = {.y = total.y - in.y, .y2 = total.y2 - in.y2};
            if (out_count > 0) {
                double in_sse = in.y2 - ((in.y * in.y) / (double)width);
                double out_sse = out.y2 - ((out.y * out.y) / (double)out_count);
                if (in_sse < 0.0) in_sse = 0.0;
                if (out_sse < 0.0) out_sse = 0.0;
                double model_sse = in_sse + out_sse;
                double explained = ref_sse - model_sse;
                if (explained < 0.0) explained = 0.0;
                if (double_is_finite_bits(explained) && explained > best_score) {
                    best_score = explained;
                    best_r2 = explained / ref_sse;
                    if (best_r2 < 0.0) best_r2 = 0.0;
                    if (best_r2 > 1.0 - 1.0e-15) best_r2 = 1.0 - 1.0e-15;
                    best_amp = fabs((in.y / (double)width) - (out.y / (double)out_count));
                    best_in_level = in.y / (double)width;
                    best_out_level = out.y / (double)out_count;
                    best_start = start;
                    best_width = width;
                }
            }

            if (start + 1 < n) {
                bls_add_sorted(buffer, sorted, start, -1.0, &in);
                bls_add_sorted(buffer, sorted, start + width, 1.0, &in);
            }
        }
    }

    if (best_start < 0) return result;

    double best_nll = lnFAP(2, 1, best_r2, n);
    float likelihood = bls_float_value(best_nll);
    if (!float_is_finite_bits(likelihood) || likelihood <= 0.0f) return result;

    if (prewhiten) {
        for (int i = 0; i < n; ++i) {
            buffer->y[i] -= (float)best_out_level;
        }
        double delta = best_in_level - best_out_level;
        for (int j = 0; j < best_width; ++j) {
            int pos = best_start + j;
            if (pos >= n) pos %= n;
            uint32_t idx = bls_unpack_index(sorted[pos]);
            buffer->y[idx] -= (float)delta;
        }
    }

    result.likelihood = likelihood;
    result.stat = bls_float_value(best_r2);
    result.amp = bls_float_value(best_amp);
    result.width = bls_float_value((double)best_width / (double)n);
    result.valid = true;
    return result;
}

typedef eval_result_t (*eval_score_fn)(buffer_t* buffer, const parameters* params, double freq, bool prewhiten);
typedef float (*eval_rank_result_fn)(const eval_result_t* result);
typedef float (*eval_rank_peak_fn)(const peak_t* peak);
typedef double (*eval_direct_step_fn)(uint32_t n, double time_span, const parameters* params);

typedef struct {
    const char* name;
    const char* stat_label;
    const char* width_label;
    const char* peak_power_label;
    eval_score_fn score;
    eval_rank_result_fn rank_result;
    eval_rank_peak_fn rank_peak;
    eval_direct_step_fn direct_frequency_step;
    bool direct_replaces_aov;
    bool uses_duration_grid;
} eval_method_t;

static inline float eval_rank_result_by_stat(const eval_result_t* result) { return result ? result->stat : 0.0f; }

static inline float eval_rank_result_by_likelihood(const eval_result_t* result) { return result ? result->likelihood : 0.0f; }

static inline float eval_rank_peak_by_stat(const peak_t* peak) { return peak ? peak->r2 : 0.0f; }

static inline float eval_rank_peak_by_likelihood(const peak_t* peak) { return peak ? peak->p : 0.0f; }

static inline double eval_gb_direct_frequency_step(uint32_t n, double time_span, const parameters* params) {
    return gb_direct_frequency_step(n, time_span, params->oversamplingFactor, params->gbAlpha);
}

static inline double eval_bls_direct_frequency_step(uint32_t n, double time_span, const parameters* params) {
    (void)n;
    return bls_direct_frequency_step(time_span, params->oversamplingFactor, params->blsMinRelWidth);
}

static inline eval_result_t eval_score_gbls(buffer_t* buffer, const parameters* params, double freq, bool prewhiten) {
    eval_result_t result = {0};
    float stat = get_r2(buffer, freq, &result.amp, prewhiten, params->gbAlpha);
    if (!float_is_finite_bits(stat)) return result;
    result.stat = stat;
    result.likelihood = (float)lnFAP(2, 1, stat, buffer->n);
    result.valid = float_is_finite_bits(result.likelihood);
    return result;
}

static inline eval_result_t eval_score_gbaw(buffer_t* buffer, const parameters* params, double freq, bool prewhiten) {
    eval_result_t result = {0};
    float stat = prewhiten ? get_r2(buffer, freq, &result.amp, true, params->gbAlpha) : get_gbaw_t(buffer, freq, &result.amp, params->gbAlpha);
    if (!float_is_finite_bits(stat)) return result;
    result.stat = stat;
    result.likelihood = prewhiten ? stat : (float)get_z(stat, buffer->n);
    result.valid = float_is_finite_bits(result.likelihood);
    return result;
}

static inline eval_result_t eval_score_bls(buffer_t* buffer, const parameters* params, double freq, bool prewhiten) {
    return get_bls_result(buffer, params, freq, prewhiten);
}

enum { EVAL_METHOD_COUNT = 3 };

static const eval_method_t EVAL_METHODS[EVAL_METHOD_COUNT] = {
    {"gbls", "R2", NULL, "NLL (Base-10)", eval_score_gbls, eval_rank_result_by_stat, eval_rank_peak_by_stat, eval_gb_direct_frequency_step, false, false},
    {"gbaw", "T", NULL, "NLL (Base-10)", eval_score_gbaw, eval_rank_result_by_stat, eval_rank_peak_by_stat, eval_gb_direct_frequency_step, false, false},
    {"bls", "R2", "Duration", "NLL (Base-10)", eval_score_bls, eval_rank_result_by_likelihood, eval_rank_peak_by_likelihood, eval_bls_direct_frequency_step,
     true, true},
};

static inline unsigned eval_method_index(gb_eval_mode mode) {
    unsigned idx = (unsigned)mode;
    return idx < EVAL_METHOD_COUNT ? idx : 0U;
}

static inline const eval_method_t* eval_method_for_mode(gb_eval_mode mode) { return &EVAL_METHODS[eval_method_index(mode)]; }

static inline const eval_method_t* eval_method_for_params(const parameters* params) { return eval_method_for_mode(params->gbEvalMode); }

static inline bool eval_method_uses_direct_grid(const eval_method_t* method, bool use_aov, int mode) {
    return mode_uses_direct_eval_grid(mode) && (!use_aov || (method && method->direct_replaces_aov));
}

static inline eval_result_t get_eval_result(buffer_t* buffer, const parameters* params, const eval_method_t* method, double freq) {
    const eval_method_t* selected = method ? method : eval_method_for_params(params);
    return selected->score(buffer, params, freq, false);
}

static inline float get_eval_likelihood(buffer_t* buffer, const parameters* params, const eval_method_t* method, double freq, float* amp, float* stat) {
    eval_result_t result = get_eval_result(buffer, params, method, freq);
    if (amp) *amp = result.amp;
    if (stat) *stat = result.stat;
    return result.valid ? result.likelihood : 0.0f;
}

static inline float eval_result_rank(const eval_method_t* method, const eval_result_t* result) {
    return method && method->rank_result ? method->rank_result(result) : 0.0f;
}

static inline float eval_peak_rank(const eval_method_t* method, const peak_t* peak) { return method && method->rank_peak ? method->rank_peak(peak) : 0.0f; }

static inline void set_peak_eval_result(peak_t* peak, const eval_result_t* result) {
    peak->p = result->likelihood;
    peak->r2 = result->stat;
    peak->amp = result->amp;
    peak->width = result->width;
}

static inline eval_result_t eval_prewhiten(buffer_t* buffer, const parameters* params, const eval_method_t* method, double freq) {
    const eval_method_t* selected = method ? method : eval_method_for_params(params);
    return selected->score(buffer, params, freq, true);
}

static inline void binsearch_peak(peak_t* peak, buffer_t* buffer, const parameters* params, const eval_method_t* method, double df) {
    const eval_method_t* selected = method ? method : eval_method_for_params(params);
    double step = df;
    eval_result_t best = get_eval_result(buffer, params, selected, peak->freq);
    if (best.valid) {
        set_peak_eval_result(peak, &best);
    } else {
        best.likelihood = peak->p;
        best.stat = peak->r2;
        best.amp = peak->amp;
    }
    float best_rank = eval_result_rank(selected, &best);

    for (int i = 0; i < 12; i++) {
        step *= 0.5;

        double freq_tmp = peak->freq - step;
        if (freq_tmp > 0.0) {
            eval_result_t candidate = get_eval_result(buffer, params, selected, freq_tmp);
            float rank = eval_result_rank(selected, &candidate);
            if (candidate.valid && rank > best_rank) {
                best_rank = rank;
                set_peak_eval_result(peak, &candidate);
                best = candidate;
                peak->freq = freq_tmp;
            }
        }

        freq_tmp = peak->freq + step;
        if (freq_tmp > 0.0) {
            eval_result_t candidate = get_eval_result(buffer, params, selected, freq_tmp);
            float rank = eval_result_rank(selected, &candidate);
            if (candidate.valid && rank > best_rank) {
                best_rank = rank;
                set_peak_eval_result(peak, &candidate);
                best = candidate;
                peak->freq = freq_tmp;
            }
        }
    }
}

static inline bool eval_peak_needs_result(const peak_t* peak) { return !float_is_finite_bits(peak->p) || !float_is_finite_bits(peak->r2) || peak->p <= 0.0f; }

static inline void evaluate_peak_at_current_frequency(peak_t* peak, buffer_t* buffer, const parameters* params, const eval_method_t* method) {
    eval_result_t result = get_eval_result(buffer, params, method, peak->freq);
    if (result.valid) {
        set_peak_eval_result(peak, &result);
    }
}

static inline double correct_ihs_res(const double sum, const int n) {
    double logSum = 1.0;
    double term = 1.0;
    for (int i = 1; i < n; i++) {
        term *= sum / (double)i;
        logSum += term;
    }
    return sum - log(logSum);
}

static inline void sortPeaks(peak_t* peaks, int length, buffer_t* buf, const parameters* params, double df) {
    const eval_method_t* method = eval_method_for_params(params);
    int i, j;
    peak_t key;

    if (params->mode == 0) {
        for (i = 0; i < length; i++) {
            peaks[i].p = correct_ihs_res(peaks[i].p, params->nterms);
        }  // Erlang's logarithmic correction'
    } else {
        for (i = 0; i < length; i++) {
            if (mode_defers_peak_evaluation(params->mode)) {
                evaluate_peak_at_current_frequency(&peaks[i], buf, params, method);
            }
            if (mode_refines_retained_peaks(params->mode)) {
                binsearch_peak(&peaks[i], buf, params, method, df);
            } else if (eval_peak_needs_result(&peaks[i])) {
                evaluate_peak_at_current_frequency(&peaks[i], buf, params, method);
            }
        };

        for (i = 1; i < length; i++) {
            key = peaks[i];
            j = i - 1;
            while (j >= 0 && eval_peak_rank(method, &peaks[j]) < eval_peak_rank(method, &key)) {
                peaks[j + 1] = peaks[j];
                j = j - 1;
            }
            peaks[j + 1] = key;
        }
    }
};

#    define FTEST_H
#endif
