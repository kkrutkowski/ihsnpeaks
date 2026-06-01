#ifndef FTEST_H

#    include <fdist.h>
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
    double L = 2 * r + 1;
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
    double norm = 1.0 / ((2 * r + 1) * (2 * r + 1) * (2 * r + 1) * (2 * r + 1));

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

    double freq_tmp = 1024.0 * freq;  // 10-bit key
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
        input[i].parts.key = ((uint16_t)(buffer->x[i] * freq_tmp)) & 0b0000001111111111;  // get rid of the 6 most significant bits
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

static inline float get_gb_stat(buffer_t* buffer, double freq, float* amp, gb_eval_mode evalMode, float gbAlpha) {
    return evalMode == GB_EVAL_GBAW ? get_gbaw_t(buffer, freq, amp, gbAlpha) : get_r2(buffer, freq, amp, false, gbAlpha);
}

static inline float get_gb_likelihood_from_stat(float stat, int n, gb_eval_mode evalMode) {
    return evalMode == GB_EVAL_GBAW ? (float)get_z(stat, n) : (float)lnFAP(2, 1, stat, n);
}

static inline float get_gb_likelihood(buffer_t* buffer, double freq, float* amp, float* stat, gb_eval_mode evalMode, float gbAlpha) {
    float value = get_gb_stat(buffer, freq, amp, evalMode, gbAlpha);
    if (stat) {
        *stat = value;
    }
    return get_gb_likelihood_from_stat(value, buffer->n, evalMode);
}

static inline void binsearch_peak(peak_t* peak, buffer_t* buffer, double df, gb_eval_mode evalMode, float gbAlpha) {
    double step = df;
    float stat;
    double freq_tmp;
    for (int i = 0; i < 12; i++) {
        step *= 0.5;

        freq_tmp = peak->freq - step;
        stat = get_gb_stat(buffer, freq_tmp, NULL, evalMode, gbAlpha);
        if (stat > peak->r2) {
            peak->r2 = stat;
            peak->freq = freq_tmp;
        }

        freq_tmp = peak->freq + step;
        stat = get_gb_stat(buffer, freq_tmp, NULL, evalMode, gbAlpha);
        if (stat > peak->r2) {
            peak->r2 = stat;
            peak->freq = freq_tmp;
        }
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

static inline void sortPeaks(peak_t* peaks, int length, buffer_t* buf, int mode, double df, int n, gb_eval_mode evalMode, float gbAlpha) {
    int i, j;
    peak_t key;

    if (mode == 0) {
        for (i = 0; i < length; i++) {
            peaks[i].p = correct_ihs_res(peaks[i].p, n);
        }  // Erlang's logarithmic correction'
    } else {
        for (i = 0; i < length; i++) {
            if (mode_refines_retained_peaks(mode)) {
                binsearch_peak(&peaks[i], buf, df, evalMode, gbAlpha);
            }
            if (mode < 4 || mode == 5) {
                peaks[i].r2 = get_gb_stat(buf, peaks[i].freq, &peaks[i].amp, evalMode, gbAlpha);
            }
            peaks[i].p = get_gb_likelihood_from_stat(peaks[i].r2, buf->n, evalMode);
        };

        for (i = 1; i < length; i++) {
            key = peaks[i];
            j = i - 1;
            while (j >= 0 && peaks[j].r2 < key.r2) {
                peaks[j + 1] = peaks[j];
                j = j - 1;
            }
            peaks[j + 1] = key;
        }
    }
};

#    define FTEST_H
#endif
