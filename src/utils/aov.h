#ifndef AOV_H
#define AOV_H

#include <fdist.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "simd.h"
#include "trig.h"

#ifndef M_SQRT1_2
#    define M_SQRT1_2 0.70710678118654752440
#endif

#define AOV_CONDITION_LIMIT 1.0e4f
#define AOV_R2_MAX (1.0f - 1.0e-7f)

typedef struct {
    float ws;
    float ymean;
    float yws;
    float chi2_ref;
} aov_reference_t;

static inline size_t aov_round_alloc(size_t size) { return (size + 63U) & ~(size_t)63U; }

static inline void *aov_aligned_alloc(size_t count, size_t size) {
    if (count == 0 || size == 0 || count > SIZE_MAX / size) return NULL;
    return aligned_alloc(64, aov_round_alloc(count * size));
}

static inline float aov_invalid_value(void) { return NAN; }

static inline bool aov_target_has_dof(const buffer_t *buffer, int degree) { return periodogram_effective_n(buffer) > (2 * degree + 1); }

static inline size_t aov_single_scratch_len(int degree) {
    if (degree <= 0) return 0U;
    return (size_t)(26 * degree + 16);
}

static inline float aov_clamp_r2(float r2) {
    if (!float_is_finite_bits(r2)) return aov_invalid_value();
    if (r2 < 0.0f) return 0.0f;
    if (r2 > AOV_R2_MAX) return AOV_R2_MAX;
    return r2;
}

static inline float aov_likelihood_from_r2(float r2, int degree, int n_eff) {
    if (n_eff <= 2 * degree + 1 || !float_is_finite_bits(r2)) return 0.0f;
    r2 = aov_clamp_r2(r2);
    if (!float_is_finite_bits(r2)) return 0.0f;
    return (float)lnFAP((2 * degree) + 1, 2 * degree, (double)r2, n_eff);
}

static inline bool aov_prepare_reference(const buffer_t *buffer, float epsilon, aov_reference_t *ref) {
    float ws = 0.0f;
    float wy = 0.0f;
    for (uint32_t i = 0; i < buffer->n; ++i) {
        float weight = 1.0f / ((buffer->dy[i] * buffer->dy[i]) + epsilon);
        ws += weight;
        wy += weight * buffer->y[i];
    }
    if (ws == 0.0f || !float_is_finite_bits(ws)) return false;

    float ymean = wy / ws;
    float yws = 0.0f;
    float chi2_ref = 0.0f;
    for (uint32_t i = 0; i < buffer->n; ++i) {
        float weight = 1.0f / ((buffer->dy[i] * buffer->dy[i]) + epsilon);
        float yc = buffer->y[i] - ymean;
        float yw = yc * weight;
        yws += yw;
        chi2_ref += yc * yw;
    }
    if (chi2_ref == 0.0f || !float_is_finite_bits(chi2_ref)) return false;

    ref->ws = ws;
    ref->ymean = ymean;
    ref->yws = yws;
    ref->chi2_ref = chi2_ref;
    return true;
}

static inline VEC aov_reflection_condition_step_vec(VEC abs2) {
    VEC out;
    for (int lane = 0; lane < VEC_LEN; ++lane) {
        float abs_reflection = sqrtf(abs2.values[lane]);
        out.values[lane] = (1.0f + abs_reflection) / (1.0f - abs_reflection);
    }
    return out;
}

static inline float aov_reflection_condition_step(float abs2) {
    float abs_reflection = sqrtf(abs2);
    return (1.0f + abs_reflection) / (1.0f - abs_reflection);
}

static inline bool aov_condition_is_singular(float condition_bound) {
    return float_is_nan_bits(condition_bound) || condition_bound < 0.0f || condition_bound > AOV_CONDITION_LIMIT;
}

static inline VEC aov_solve_levinson_vec(size_t n, const VEC *restrict Rr, const VEC *restrict Ri, const VEC *restrict Yr, const VEC *restrict Yi,
                                         VEC *restrict Xr, VEC *restrict Xi, VEC *restrict Ar, VEC *restrict Ai, VEC *restrict Apr, VEC *restrict Api) {
    VEC E = Rr[0];
    VEC cond_bound = SET_VEC(1.0f);

    Xr[0].data = Yr[0].data / E.data;
    Xi[0].data = Yi[0].data / E.data;
    Ar[0] = SET_VEC(1.0f);
    Ai[0] = SET_VEC(0.0f);

    for (size_t k = 1; k < n; ++k) {
        VEC lambda_r = SET_VEC(0.0f);
        VEC lambda_i = SET_VEC(0.0f);
        for (size_t i = 0; i < k; ++i) {
            lambda_r.data += Rr[k - i].data * Ar[i].data - Ri[k - i].data * Ai[i].data;
            lambda_i.data += Rr[k - i].data * Ai[i].data + Ri[k - i].data * Ar[i].data;
        }

        VEC gamma_r = {.data = -lambda_r.data / E.data};
        VEC gamma_i = {.data = -lambda_i.data / E.data};

        for (size_t i = 0; i < k; ++i) {
            Apr[i] = Ar[i];
            Api[i] = Ai[i];
        }

        Ar[k].data = gamma_r.data * Apr[0].data + gamma_i.data * Api[0].data;
        Ai[k].data = -gamma_r.data * Api[0].data + gamma_i.data * Apr[0].data;

        for (size_t i = 1; i < k; ++i) {
            VEC ap_r = Apr[k - i];
            VEC ap_i = {.data = -Api[k - i].data};
            Ar[i].data = Apr[i].data + gamma_r.data * ap_r.data - gamma_i.data * ap_i.data;
            Ai[i].data = Api[i].data + gamma_r.data * ap_i.data + gamma_i.data * ap_r.data;
        }

        VEC abs2_gamma = {.data = gamma_r.data * gamma_r.data + gamma_i.data * gamma_i.data};
        cond_bound.data *= aov_reflection_condition_step_vec(abs2_gamma).data;
        E.data *= SET_VEC(1.0f).data - abs2_gamma.data;

        VEC mu_r = Yr[k];
        VEC mu_i = Yi[k];
        for (size_t i = 0; i < k; ++i) {
            mu_r.data -= Rr[k - i].data * Xr[i].data - Ri[k - i].data * Xi[i].data;
            mu_i.data -= Rr[k - i].data * Xi[i].data + Ri[k - i].data * Xr[i].data;
        }

        VEC nu_r = {.data = mu_r.data / E.data};
        VEC nu_i = {.data = mu_i.data / E.data};

        Xr[k] = nu_r;
        Xi[k] = nu_i;
        for (size_t i = 0; i < k; ++i) {
            VEC ak_r = Ar[k - i];
            VEC ak_i = {.data = -Ai[k - i].data};
            Xr[i].data += nu_r.data * ak_r.data - nu_i.data * ak_i.data;
            Xi[i].data += nu_r.data * ak_i.data + nu_i.data * ak_r.data;
        }
    }

    return cond_bound;
}

static inline float aov_solve_levinson_scalar(size_t n, const float *restrict Rr, const float *restrict Ri, const float *restrict Yr, const float *restrict Yi,
                                              float *restrict Xr, float *restrict Xi, float *restrict Ar, float *restrict Ai, float *restrict Apr,
                                              float *restrict Api) {
    float E = Rr[0];
    float cond_bound = 1.0f;
    Xr[0] = Yr[0] / E;
    Xi[0] = Yi[0] / E;
    Ar[0] = 1.0f;
    Ai[0] = 0.0f;

    for (size_t k = 1; k < n; ++k) {
        float lambda_r = 0.0f;
        float lambda_i = 0.0f;
        for (size_t i = 0; i < k; ++i) {
            lambda_r += Rr[k - i] * Ar[i] - Ri[k - i] * Ai[i];
            lambda_i += Rr[k - i] * Ai[i] + Ri[k - i] * Ar[i];
        }

        float gamma_r = -lambda_r / E;
        float gamma_i = -lambda_i / E;

        for (size_t i = 0; i < k; ++i) {
            Apr[i] = Ar[i];
            Api[i] = Ai[i];
        }

        Ar[k] = gamma_r * Apr[0] + gamma_i * Api[0];
        Ai[k] = -gamma_r * Api[0] + gamma_i * Apr[0];

        for (size_t i = 1; i < k; ++i) {
            float ap_r = Apr[k - i];
            float ap_i = -Api[k - i];
            Ar[i] = Apr[i] + gamma_r * ap_r - gamma_i * ap_i;
            Ai[i] = Api[i] + gamma_r * ap_i + gamma_i * ap_r;
        }

        float abs2_gamma = gamma_r * gamma_r + gamma_i * gamma_i;
        cond_bound *= aov_reflection_condition_step(abs2_gamma);
        E *= 1.0f - abs2_gamma;

        float mu_r = Yr[k];
        float mu_i = Yi[k];
        for (size_t i = 0; i < k; ++i) {
            mu_r -= Rr[k - i] * Xr[i] - Ri[k - i] * Xi[i];
            mu_i -= Rr[k - i] * Xi[i] + Ri[k - i] * Xr[i];
        }

        float nu_r = mu_r / E;
        float nu_i = mu_i / E;
        Xr[k] = nu_r;
        Xi[k] = nu_i;
        for (size_t i = 0; i < k; ++i) {
            float ak_r = Ar[k - i];
            float ak_i = -Ai[k - i];
            Xr[i] += nu_r * ak_r - nu_i * ak_i;
            Xi[i] += nu_r * ak_i + nu_i * ak_r;
        }
    }

    return cond_bound;
}

static inline void aov_gls_impl(const float *Sw, const float *Cw, const float *Syw, const float *Cyw, int nfreq, float ws, float yws, float chi2_ref,
                                float *power) {
    float inv_ws = 1.0f / ws;
    float inv_chi2_ref = 1.0f / chi2_ref;
    float half = 0.5f;
    float sqrt_half = (float)M_SQRT1_2;

    for (int idx = 0; idx < nfreq; ++idx) {
        float S = Sw[(size_t)nfreq + idx];
        float C = Cw[(size_t)nfreq + idx];
        float S2 = Sw[(size_t)2 * (size_t)nfreq + idx];
        float C2 = Cw[(size_t)2 * (size_t)nfreq + idx];
        float Sh = Syw[(size_t)nfreq + idx];
        float Ch = Cyw[(size_t)nfreq + idx];

        float num = S2 - (2.0f * S * C * inv_ws);
        float den = C2 - ((C * C - S * S) * inv_ws);
        float radius = sqrtf((num * num) + (den * den));

        float S2w = 0.0f;
        float C2w = 1.0f;
        if (radius != 0.0f) {
            float den_abs = den;
            float den_sign = 1.0f;
            if (den < 0.0f) {
                den_abs = -den;
                den_sign = -1.0f;
            }
            S2w = den_sign * num / radius;
            C2w = den_abs / radius;
        }

        float one_plus_c2 = 1.0f + C2w;
        float one_minus_c2 = 1.0f - C2w;
        if (one_plus_c2 < 0.0f) one_plus_c2 = 0.0f;
        if (one_minus_c2 < 0.0f) one_minus_c2 = 0.0f;

        float Ctau = sqrt_half * sqrtf(one_plus_c2);
        float Stau = 0.0f;
        if (S2w > 0.0f)
            Stau = sqrt_half * sqrtf(one_minus_c2);
        else if (S2w < 0.0f)
            Stau = -sqrt_half * sqrtf(one_minus_c2);

        float YC = (Ch * Ctau) + (Sh * Stau);
        float YS = (Sh * Ctau) - (Ch * Stau);
        float XC = (C * Ctau) + (S * Stau);
        float XS = (S * Ctau) - (C * Stau);

        float CC = half * (ws + ((C2 * C2w) + (S2 * S2w)));
        float SS = half * (ws - ((C2 * C2w) + (S2 * S2w)));
        CC -= (XC * XC) * inv_ws;
        SS -= (XS * XS) * inv_ws;
        YC -= (yws * XC) * inv_ws;
        YS -= (yws * XS) * inv_ws;

        power[idx] = aov_clamp_r2((((YC * YC) / CC) + ((YS * YS) / SS)) * inv_chi2_ref);
    }
}

static inline float aov_solve_single_from_sums(const float *Sw, const float *Cw, const float *Syw, const float *Cyw, int degree, float chi2_ref, float *scratch,
                                               size_t scratch_len) {
    if (degree == 1) {
        float power = aov_invalid_value();
        aov_gls_impl(Sw, Cw, Syw, Cyw, 1, Cw[0], Cyw[0], chi2_ref, &power);
        return power;
    }

    int norder = (2 * degree) + 1;
    size_t required = ((size_t)norder + 1U) * 2U + ((size_t)norder * 8U);
    if (!scratch || scratch_len < required) return aov_invalid_value();

    size_t offset = 0;
    float *Rr = scratch + offset;
    offset += (size_t)norder + 1U;
    float *Ri = scratch + offset;
    offset += (size_t)norder + 1U;
    float *Yr = scratch + offset;
    offset += (size_t)norder;
    float *Yi = scratch + offset;
    offset += (size_t)norder;
    float *Xr = scratch + offset;
    offset += (size_t)norder;
    float *Xi = scratch + offset;
    offset += (size_t)norder;
    float *Ar = scratch + offset;
    offset += (size_t)norder;
    float *Ai = scratch + offset;
    offset += (size_t)norder;
    float *Apr = scratch + offset;
    offset += (size_t)norder;
    float *Api = scratch + offset;

    for (int k = 0; k < norder; ++k) {
        Rr[k] = Cw[k];
        Ri[k] = -Sw[k];
        int h = k - degree;
        if (h == 0) {
            Yr[k] = Cyw[0];
            Yi[k] = 0.0f;
        } else if (h > 0) {
            Yr[k] = Cyw[h];
            Yi[k] = -Syw[h];
        } else {
            Yr[k] = Cyw[-h];
            Yi[k] = Syw[-h];
        }
    }
    Rr[norder] = 0.0f;
    Ri[norder] = 0.0f;

    float condition = aov_solve_levinson_scalar((size_t)norder, Rr, Ri, Yr, Yi, Xr, Xi, Ar, Ai, Apr, Api);
    if (aov_condition_is_singular(condition)) return aov_invalid_value();

    float dot = 0.0f;
    for (int k = 0; k < norder; ++k) dot += (Yr[k] * Xr[k]) + (Yi[k] * Xi[k]);
    return aov_clamp_r2(dot / chi2_ref);
}

static inline int aov_solve_periodogram_vec(const float *Sw, const float *Cw, const float *Syw, const float *Cyw, int nfreq, int degree, float chi2_ref,
                                            float *power) {
    int norder = (2 * degree) + 1;
    VEC *Rr = (VEC *)aov_aligned_alloc((size_t)norder + 1U, sizeof(VEC));
    VEC *Ri = (VEC *)aov_aligned_alloc((size_t)norder + 1U, sizeof(VEC));
    VEC *Yr = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Yi = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Xr = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Xi = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Ar = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Ai = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Apr = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    VEC *Api = (VEC *)aov_aligned_alloc((size_t)norder, sizeof(VEC));
    if (!Rr || !Ri || !Yr || !Yi || !Xr || !Xi || !Ar || !Ai || !Apr || !Api) {
        free(Rr);
        free(Ri);
        free(Yr);
        free(Yi);
        free(Xr);
        free(Xi);
        free(Ar);
        free(Ai);
        free(Apr);
        free(Api);
        return -1;
    }

    float inv_chi2_ref = 1.0f / chi2_ref;
    for (int base = 0; base < nfreq; base += VEC_LEN) {
        for (int k = 0; k < norder; ++k) {
            int h = k - degree;
            for (int lane = 0; lane < VEC_LEN; ++lane) {
                int idx = base + lane;
                if (idx >= nfreq) idx = nfreq - 1;
                Rr[k].values[lane] = Cw[(size_t)k * (size_t)nfreq + (size_t)idx];
                Ri[k].values[lane] = -Sw[(size_t)k * (size_t)nfreq + (size_t)idx];
                if (h == 0) {
                    Yr[k].values[lane] = Cyw[idx];
                    Yi[k].values[lane] = 0.0f;
                } else if (h > 0) {
                    Yr[k].values[lane] = Cyw[(size_t)h * (size_t)nfreq + (size_t)idx];
                    Yi[k].values[lane] = -Syw[(size_t)h * (size_t)nfreq + (size_t)idx];
                } else {
                    int ah = -h;
                    Yr[k].values[lane] = Cyw[(size_t)ah * (size_t)nfreq + (size_t)idx];
                    Yi[k].values[lane] = Syw[(size_t)ah * (size_t)nfreq + (size_t)idx];
                }
            }
        }
        Rr[norder] = SET_VEC(0.0f);
        Ri[norder] = SET_VEC(0.0f);

        VEC condition = aov_solve_levinson_vec((size_t)norder, Rr, Ri, Yr, Yi, Xr, Xi, Ar, Ai, Apr, Api);
        for (int lane = 0; lane < VEC_LEN; ++lane) {
            int idx = base + lane;
            if (idx >= nfreq) continue;
            if (aov_condition_is_singular(condition.values[lane])) {
                power[idx] = aov_invalid_value();
                continue;
            }
            float dot = 0.0f;
            for (int k = 0; k < norder; ++k) dot += (Yr[k].values[lane] * Xr[k].values[lane]) + (Yi[k].values[lane] * Xi[k].values[lane]);
            power[idx] = aov_clamp_r2(dot * inv_chi2_ref);
        }
    }

    free(Rr);
    free(Ri);
    free(Yr);
    free(Yi);
    free(Xr);
    free(Xi);
    free(Ar);
    free(Ai);
    free(Apr);
    free(Api);
    return 0;
}

static inline void aov_fill_complex_input(buffer_t *buffer, double base_freq, double freq_factor, float epsilon, float ymean, bool yweighted) {
    uint32_t i = 0;
    for (; i < buffer->n; ++i) {
        float weight = 1.0f / ((buffer->dy[i] * buffer->dy[i]) + epsilon);
        float value = yweighted ? weight * (buffer->y[i] - ymean) : weight;
        float phase = (float)(buffer->x[i] * base_freq * freq_factor);
        buffer->inputReal[i] = value * cos2pif_tls(phase);
        buffer->inputImag[i] = value * sin2pif_tls(phase);
    }
    for (; i < buffer->paddedLen; ++i) {
        buffer->inputReal[i] = 0.0f;
        buffer->inputImag[i] = 0.0f;
    }
}

static inline void aov_compute_trig_sums_block(buffer_t *buffer, parameters *params, double block_fmin, uint32_t count, int max_factor, float epsilon,
                                               float ymean, bool yweighted, float *S, float *C) {
    for (int q = 1; q <= max_factor; ++q) {
        double freq_factor = (double)q;
        aov_fill_complex_input(buffer, block_fmin, freq_factor, epsilon, ymean, yweighted);
        nufft1_execute(buffer->nufftWorkspace, buffer->inputReal, buffer->inputImag, buffer->blockReal, buffer->blockImag, q);
        memcpy(C + ((size_t)q * (size_t)count), buffer->blockReal, (size_t)count * sizeof(float));
        memcpy(S + ((size_t)q * (size_t)count), buffer->blockImag, (size_t)count * sizeof(float));
    }
}

static inline bool aov_quadratic_peak_position(double freq, double fstep, float left, float center, float right, float oversampling, double *peak_freq,
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

typedef struct {
    bool has_left;
    bool has_center;
    uint32_t left_idx;
    uint32_t center_idx;
    double left_freq;
    double center_freq;
    float left_value;
    float center_value;
} aov_peak_stream_t;

static inline void aov_append_peak(buffer_t *buffer, const parameters *params, double freq, float r2, double df);

static inline void aov_stream_peak_value(aov_peak_stream_t *stream, buffer_t *buffer, const parameters *params, double freq, float value, uint32_t idx,
                                         uint32_t nfreq, float threshold, double fstep, double df, char *stringBuff, int precision, bool write_spectrum) {
    if (write_spectrum) {
        float magnitude = aov_likelihood_from_r2(value, params->nterms, periodogram_effective_n(buffer));
        if (!float_is_finite_bits(magnitude)) magnitude = 0.0f;
        appendFreq(freq, magnitude, precision, params->outputPeriod, &buffer->spectrum, stringBuff);
    }

    if (stream->has_left && stream->has_center && stream->center_idx > 0U && idx < nfreq) {
        double peak_freq = stream->center_freq;
        float peak_magnitude = stream->center_value;
        if (aov_quadratic_peak_position(stream->center_freq, fstep, stream->left_value, stream->center_value, value, params->oversamplingFactor, &peak_freq,
                                        &peak_magnitude) &&
            stream->center_value > stream->left_value && stream->center_value > value &&
            (stream->center_value > threshold || mode_evaluates_all_local_peaks(params->mode))) {
            aov_append_peak(buffer, params, peak_freq, peak_magnitude, df);
        }
    }

    stream->has_left = stream->has_center;
    stream->left_idx = stream->center_idx;
    stream->left_freq = stream->center_freq;
    stream->left_value = stream->center_value;
    stream->has_center = true;
    stream->center_idx = idx;
    stream->center_freq = freq;
    stream->center_value = value;
}

static inline int execute_aov_sweep(buffer_t *buffer, parameters *params, double fmin, double fstep, uint32_t nfreq, float threshold, double df,
                                    bool write_spectrum, int precision, char *stringBuff) {
    if (nfreq == 0) return 0;
    if (!aov_target_has_dof(buffer, params->nterms)) {
        if (write_spectrum) {
            for (uint32_t i = 0; i < nfreq; ++i) {
                double freq = fmin + ((double)i * fstep);
                appendFreq(freq, 0.0f, precision, params->outputPeriod, &buffer->spectrum, stringBuff);
            }
        }
        return 0;
    }

    aov_reference_t ref;
    if (!aov_prepare_reference(buffer, params->epsilon, &ref)) {
        if (write_spectrum) {
            for (uint32_t i = 0; i < nfreq; ++i) {
                double freq = fmin + ((double)i * fstep);
                appendFreq(freq, 0.0f, precision, params->outputPeriod, &buffer->spectrum, stringBuff);
            }
        }
        return 0;
    }

    int degree = params->nterms;
    int max_factor = 2 * degree;
    uint32_t block_len = buffer->activeOutputLen;
    float *Sw = (float *)aov_aligned_alloc((size_t)(max_factor + 1) * (size_t)block_len, sizeof(float));
    float *Cw = (float *)aov_aligned_alloc((size_t)(max_factor + 1) * (size_t)block_len, sizeof(float));
    float *Syw = (float *)aov_aligned_alloc((size_t)(degree + 1) * (size_t)block_len, sizeof(float));
    float *Cyw = (float *)aov_aligned_alloc((size_t)(degree + 1) * (size_t)block_len, sizeof(float));
    float *block_power = (float *)aov_aligned_alloc((size_t)block_len, sizeof(float));
    if (!Sw || !Cw || !Syw || !Cyw || !block_power) {
        free(Sw);
        free(Cw);
        free(Syw);
        free(Cyw);
        free(block_power);
        return -1;
    }

    nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);

    size_t num_blocks = ((size_t)nfreq + (size_t)block_len - 1U) / (size_t)block_len;
    int status = 0;
    aov_peak_stream_t stream = {0};

    for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
        uint32_t base = (uint32_t)(block_idx * (size_t)block_len);
        uint32_t count = block_len;
        if (base + count > nfreq) count = nfreq - base;
        double block_fmin = fmin + ((double)base * fstep);

        for (uint32_t i = 0; i < count; ++i) {
            Sw[i] = 0.0f;
            Cw[i] = ref.ws;
            Syw[i] = 0.0f;
            Cyw[i] = ref.yws;
        }

        aov_compute_trig_sums_block(buffer, params, block_fmin, count, max_factor, params->epsilon, ref.ymean, false, Sw, Cw);
        aov_compute_trig_sums_block(buffer, params, block_fmin, count, degree, params->epsilon, ref.ymean, true, Syw, Cyw);

        if (degree == 1) {
            aov_gls_impl(Sw, Cw, Syw, Cyw, (int)count, ref.ws, ref.yws, ref.chi2_ref, block_power);
        } else {
            status = aov_solve_periodogram_vec(Sw, Cw, Syw, Cyw, (int)count, degree, ref.chi2_ref, block_power);
            if (status != 0) break;
        }

        for (uint32_t k = 0; k < count; ++k) {
            uint32_t idx = base + k;
            double freq = fmin + ((double)idx * fstep);
            if (write_spectrum) buffer->power[idx] = block_power[k];
            aov_stream_peak_value(&stream, buffer, params, freq, block_power[k], idx, nfreq, threshold, fstep, df, stringBuff, precision, write_spectrum);
        }
    }

    free(Sw);
    free(Cw);
    free(Syw);
    free(Cyw);
    free(block_power);
    return status;
}

static inline float aov_get_stat(buffer_t *buffer, const parameters *params, double freq) {
    int degree = params->nterms;
    if (!aov_target_has_dof(buffer, degree)) return aov_invalid_value();

    aov_reference_t ref;
    if (!aov_prepare_reference(buffer, params->epsilon, &ref)) return aov_invalid_value();

    size_t scratch_len = aov_single_scratch_len(degree);
    if (!buffer->aovScratch || buffer->aovScratchLen < scratch_len) return aov_invalid_value();
    float *scratch = buffer->aovScratch;
    size_t offset = 0;
    int trig_len = (2 * degree) + 1;
    int ytrig_len = degree + 1;
    float *Sw = scratch + offset;
    offset += (size_t)trig_len;
    float *Cw = scratch + offset;
    offset += (size_t)trig_len;
    float *Syw = scratch + offset;
    offset += (size_t)ytrig_len;
    float *Cyw = scratch + offset;
    offset += (size_t)ytrig_len;

    Cw[0] = ref.ws;
    Sw[0] = 0.0f;
    Cyw[0] = ref.yws;
    Syw[0] = 0.0f;

    for (int k = 1; k <= 2 * degree; ++k) {
        Cw[k] = 0.0f;
        Sw[k] = 0.0f;
    }
    for (int k = 1; k <= degree; ++k) {
        Cyw[k] = 0.0f;
        Syw[k] = 0.0f;
    }

    for (uint32_t i = 0; i < buffer->n; ++i) {
        float weight = 1.0f / ((buffer->dy[i] * buffer->dy[i]) + params->epsilon);
        float yw = weight * (buffer->y[i] - ref.ymean);
        float phase = (float)(buffer->x[i] * freq);
        for (int k = 1; k <= 2 * degree; ++k) {
            float kphase = phase * (float)k;
            float c = cos2pif_tls(kphase);
            float s = sin2pif_tls(kphase);
            Cw[k] += weight * c;
            Sw[k] += weight * s;
            if (k <= degree) {
                Cyw[k] += yw * c;
                Syw[k] += yw * s;
            }
        }
    }

    return aov_solve_single_from_sums(Sw, Cw, Syw, Cyw, degree, ref.chi2_ref, scratch + offset, buffer->aovScratchLen - offset);
}

static inline void aov_binsearch_peak(peak_t *peak, buffer_t *buffer, const parameters *params, double df) {
    float center = aov_get_stat(buffer, params, peak->freq);
    if (float_is_finite_bits(center)) peak->r2 = center;
    double step = df;
    for (int i = 0; i < 12; ++i) {
        step *= 0.5;
        double freq_tmp = peak->freq - step;
        float stat = aov_get_stat(buffer, params, freq_tmp);
        if (float_is_finite_bits(stat) && stat > peak->r2) {
            peak->r2 = stat;
            peak->freq = freq_tmp;
        }
        freq_tmp = peak->freq + step;
        stat = aov_get_stat(buffer, params, freq_tmp);
        if (float_is_finite_bits(stat) && stat > peak->r2) {
            peak->r2 = stat;
            peak->freq = freq_tmp;
        }
    }
}

static inline void aov_append_peak(buffer_t *buffer, const parameters *params, double freq, float r2, double df) {
    if (params->npeaks <= 0) return;
    if (!float_is_finite_bits(r2)) return;
    const eval_method_t *method = eval_method_for_params(params);
    peak_t appended = {0};
    appended.freq = freq;
    appended.r2 = r2;
    appended.p = r2;

    float rank = r2;
    if (!mode_defers_peak_evaluation(params->mode)) {
        if (mode_eagerly_refines_peaks(params->mode)) {
            binsearch_peak(&appended, buffer, params, method, df);
        } else {
            evaluate_peak_at_current_frequency(&appended, buffer, params, method);
        }
        rank = eval_peak_rank(method, &appended);
        if (!float_is_finite_bits(rank)) return;
    }

    int idx = (int)buffer->nPeaks;
    float last_rank =
        mode_defers_peak_evaluation(params->mode) ? buffer->peaks[params->npeaks - 1].p : eval_peak_rank(method, &buffer->peaks[params->npeaks - 1]);
    if (idx >= params->npeaks && rank <= last_rank) return;
    while (idx > 0) {
        float prev_rank = mode_defers_peak_evaluation(params->mode) ? buffer->peaks[idx - 1].p : eval_peak_rank(method, &buffer->peaks[idx - 1]);
        if (rank <= prev_rank) break;
        --idx;
    }
    if (buffer->nPeaks < (uint32_t)params->npeaks) ++buffer->nPeaks;
    for (int i = (int)buffer->nPeaks - 1; i > idx; --i) buffer->peaks[i] = buffer->peaks[i - 1];
    if (idx < params->npeaks) buffer->peaks[idx] = appended;
}

static inline void aov_sort_peaks(peak_t *peaks, int length, buffer_t *buffer, const parameters *params, double df) {
    const eval_method_t *method = eval_method_for_params(params);
    int n_eff = periodogram_effective_n(buffer);
    if (params->mode == 0) {
        for (int i = 0; i < length; ++i) peaks[i].p = aov_likelihood_from_r2(peaks[i].r2, params->nterms, n_eff);
        return;
    }

    for (int i = 0; i < length; ++i) {
        if (mode_defers_peak_evaluation(params->mode)) evaluate_peak_at_current_frequency(&peaks[i], buffer, params, method);
        if (mode_refines_retained_peaks(params->mode)) {
            binsearch_peak(&peaks[i], buffer, params, method, df);
        } else if (eval_peak_needs_result(&peaks[i])) {
            evaluate_peak_at_current_frequency(&peaks[i], buffer, params, method);
        }
    }

    for (int i = 1; i < length; ++i) {
        peak_t key = peaks[i];
        int j = i - 1;
        while (j >= 0 && eval_peak_rank(method, &peaks[j]) < eval_peak_rank(method, &key)) {
            peaks[j + 1] = peaks[j];
            --j;
        }
        peaks[j + 1] = key;
    }
}

#endif
