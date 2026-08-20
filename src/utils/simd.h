#ifndef SIMD_H
#define SIMD_H

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#ifndef VEC_BYTES
#    ifdef __AVX512F__
#        define VEC_BYTES 64
#    elif defined(__AVX__)
#        define VEC_BYTES 32
#    else
#        define VEC_BYTES 16
#    endif
#endif

#define VEC_LEN (VEC_BYTES / (int)sizeof(float))
#define DVEC_LEN (VEC_BYTES / (int)sizeof(double))
#define IVEC_LEN (VEC_BYTES / (int)sizeof(int32_t))
#define LVEC_LEN (VEC_BYTES / (int)sizeof(int64_t))

typedef float vecf_data __attribute__((vector_size(VEC_BYTES)));
typedef double vecd_data __attribute__((vector_size(VEC_BYTES)));
typedef int32_t veci_data __attribute__((vector_size(VEC_BYTES)));
typedef int64_t vecl_data __attribute__((vector_size(VEC_BYTES)));

typedef union {
    vecf_data data;
    vecf_data gdata;
    float values[VEC_LEN];
} VEC;

typedef union {
    vecd_data data;
    vecd_data gdata;
    double values[DVEC_LEN];
} DVEC;

typedef union {
    veci_data data;
    veci_data gdata;
    int32_t values[IVEC_LEN];
} IVEC;

typedef union {
    vecl_data data;
    vecl_data gdata;
    int64_t values[LVEC_LEN];
} LVEC;

#define SET_VEC(val) ((VEC){.data = (vecf_data){} + (float)(val)})
#define SET_DVEC(val) ((DVEC){.data = (vecd_data){} + (double)(val)})
#define SET_IVEC(val) ((IVEC){.data = (veci_data){} + (int32_t)(val)})
#define SET_LVEC(val) ((LVEC){.data = (vecl_data){} + (int64_t)(val)})

#ifndef AOV_R2_MAX
#    define AOV_R2_MAX (1.0f - 1.0e-7f)
#endif

#ifndef float_is_finite_bits
static inline bool float_is_finite_bits_simd(float v) {
    union {
        float f;
        uint32_t u;
    } cast = {.f = v};
    return (cast.u & 0x7F800000U) != 0x7F800000U;
}
#    define float_is_finite_bits(v) float_is_finite_bits_simd(v)
#endif

#ifndef double_is_finite_bits
static inline bool double_is_finite_bits_simd(double v) {
    union {
        double d;
        uint64_t u;
    } cast = {.d = v};
    return (cast.u & 0x7FF0000000000000ULL) != 0x7FF0000000000000ULL;
}
#    define double_is_finite_bits(v) double_is_finite_bits_simd(v)
#endif

/* ========================================================================= */
/* Fast SIMD Logarithm Implementations (Float32 and Float64)                 */
/* ========================================================================= */

static inline VEC ln_ps(const VEC x) {
    const IVEC mant_mask = SET_IVEC(0x807FFFFF);
    const VEC ln2_vec = SET_VEC(0.69314718f);
    const IVEC exp_bias = SET_IVEC(127);
    const IVEC set_1 = SET_IVEC(0x3F800000);
    const VEC c[6] = {SET_VEC(-1.936759742e0f), SET_VEC(3.514087297e0f),  SET_VEC(-2.440029763e0f),
                      SET_VEC(1.116090027e0f),  SET_VEC(-2.83826848e-1f), SET_VEC(3.04490045e-2f)};

    union {
        vecf_data f;
        veci_data i;
    } x_cast = {.f = x.gdata};
    IVEC x_bits = {.data = x_cast.i};
    IVEC exp_bits = {.data = (x_bits.data >> 23) & SET_IVEC(0xFF).data};
    IVEC unbiased_exp = {.data = exp_bits.data - exp_bias.data};
    IVEC mant_bits = {.data = (x_bits.data & mant_mask.data) | set_1.data};
    union {
        veci_data i;
        vecf_data f;
    } mant_cast = {.i = mant_bits.gdata};
    VEC mant_vec = {.data = mant_cast.f};

    VEC exp_float = {.data = __builtin_convertvector(unbiased_exp.data, vecf_data)};
    VEC ln_mant = c[5];
    ln_mant.data = ln_mant.data * mant_vec.data + c[4].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[3].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[2].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[1].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[0].data;

    return (VEC){.data = exp_float.data * ln2_vec.data + ln_mant.data};
}

static inline DVEC ln_pd(const DVEC x) {
    const LVEC mant_mask = SET_LVEC(0x000FFFFFFFFFFFFFULL);
    const DVEC ln2_vec = SET_DVEC(0.693147180559945309417232);
    const LVEC exp_bias = SET_LVEC(1023LL);
    const LVEC set_1 = SET_LVEC(0x3FF0000000000000LL);
    const DVEC c[6] = {SET_DVEC(-1.936759742), SET_DVEC(3.514087297),  SET_DVEC(-2.440029763),
                       SET_DVEC(1.116090027),  SET_DVEC(-0.283826848), SET_DVEC(0.0304490045)};

    union {
        vecd_data d;
        vecl_data l;
    } x_cast = {.d = x.gdata};
    LVEC x_bits = {.data = x_cast.l};
    LVEC exp_bits = {.data = (x_bits.data >> 52) & SET_LVEC(0x7FFLL).data};
    LVEC unbiased_exp = {.data = exp_bits.data - exp_bias.data};
    LVEC mant_bits = {.data = (x_bits.data & mant_mask.data) | set_1.data};
    union {
        vecl_data l;
        vecd_data d;
    } mant_cast = {.l = mant_bits.data};
    DVEC mant_vec = {.data = mant_cast.d};

    DVEC exp_double = {.data = __builtin_convertvector(unbiased_exp.data, vecd_data)};
    DVEC ln_mant = c[5];
    ln_mant.data = ln_mant.data * mant_vec.data + c[4].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[3].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[2].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[1].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[0].data;

    return (DVEC){.data = exp_double.data * ln2_vec.data + ln_mant.data};
}

/* ========================================================================= */
/* IHS Rayleigh Power Correction                                             */
/* ========================================================================= */

static inline VEC correctPower(const VEC K, const float nInv) {
    const VEC n = SET_VEC(nInv);
    VEC term1 = {.data = ((SET_VEC(2.0f).data * K.data) - (K.data * K.data)) * (SET_VEC(0.25f).data * n.data)};
    VEC term2 = {.data = ((SET_VEC(24.0f).data * K.data) - (SET_VEC(132.0f).data * K.data * K.data) + (SET_VEC(76.0f).data * K.data * K.data * K.data) -
                          (SET_VEC(9.0f).data * K.data * K.data * K.data * K.data)) *
                         (n.data * n.data * SET_VEC(3.4722222e-3f).data)};
    VEC inside_log = {.data = SET_VEC(1.0f).data + term1.data - term2.data};
    VEC log_result = ln_ps(inside_log);
    return (VEC){.data = K.data - log_result.data};
}

/* ========================================================================= */
/* Robust Vectorized Multi-Harmonic Pochhammer Loops (Float32 and Float64)   */
/* ========================================================================= */

static inline VEC nll_exact_pochhammer_vec(int d, const VEC R2, float b) {
    VEC out;
    for (int lane = 0; lane < VEC_LEN; ++lane) {
        float r2 = R2.values[lane];
        if (!float_is_finite_bits(r2) || r2 <= 0.0f) {
            out.values[lane] = 0.0f;
            continue;
        }
        if (r2 > AOV_R2_MAX) r2 = AOV_R2_MAX;

        float term0;
        if (r2 < 0.15f) {
            /* Analytical polynomial for -b * ln1p(-r2) avoids precision loss near zero */
            float log_omr2 = r2 * (1.0f + r2 * (0.5f + r2 * (0.33333334f + r2 * (0.25f + r2 * (0.2f + r2 * 0.16666667f)))));
            term0 = b * log_omr2;
        } else {
            term0 = -b * logf(1.0f - r2);
        }

        if (d <= 1) {
            out.values[lane] = (term0 > 0.0f && float_is_finite_bits(term0)) ? term0 : 0.0f;
            continue;
        }

        float c = 1.0f;
        float s_minus_1 = 0.0f;
        for (int k = 1; k < d; ++k) {
            c *= (b + (float)(k - 1)) / (float)k * r2;
            s_minus_1 += c;
        }

        float ln_s;
        if (s_minus_1 < 0.15f) {
            float u = s_minus_1;
            ln_s = u * (1.0f - u * (0.5f - u * (0.33333334f - u * 0.25f)));
        } else {
            ln_s = logf(1.0f + s_minus_1);
        }

        float nll = term0 - ln_s;
        out.values[lane] = (nll > 0.0f && float_is_finite_bits(nll)) ? nll : 0.0f;
    }
    return out;
}

static inline DVEC nll_exact_pochhammer_dvec(int d, const DVEC R2, double b) {
    DVEC out;
    for (int lane = 0; lane < DVEC_LEN; ++lane) {
        double r2 = R2.values[lane];
        if (!double_is_finite_bits(r2) || r2 <= 0.0) {
            out.values[lane] = 0.0;
            continue;
        }
        if (r2 > (double)AOV_R2_MAX) r2 = (double)AOV_R2_MAX;

        double term0;
        if (r2 < 0.15) {
            double log_omr2 = r2 * (1.0 + r2 * (0.5 + r2 * (1.0 / 3.0 + r2 * (0.25 + r2 * (0.2 + r2 / 6.0)))));
            term0 = b * log_omr2;
        } else {
            term0 = -b * log(1.0 - r2);
        }

        if (d <= 1) {
            out.values[lane] = (term0 > 0.0 && double_is_finite_bits(term0)) ? term0 : 0.0;
            continue;
        }

        double c = 1.0;
        double s_minus_1 = 0.0;
        for (int k = 1; k < d; ++k) {
            c *= (b + (double)(k - 1)) / (double)k * r2;
            s_minus_1 += c;
        }

        double ln_s;
        if (s_minus_1 < 0.15) {
            double u = s_minus_1;
            ln_s = u * (1.0 - u * (0.5 - u * (1.0 / 3.0 - u * 0.25)));
        } else {
            ln_s = log(1.0 + s_minus_1);
        }

        double nll = term0 - ln_s;
        out.values[lane] = (nll > 0.0 && double_is_finite_bits(nll)) ? nll : 0.0;
    }
    return out;
}

/* ========================================================================= */
/* Robust Vectorized GBLS Continued Fractions (Float32 and Float64)          */
/* ========================================================================= */

static inline float fast_lbetaf_simd(float a, float b) { return lgammaf(a) + lgammaf(b) - lgammaf(a + b); }

static inline double fast_lbeta_simd(double a, double b) { return lgamma(a) + lgamma(b) - lgamma(a + b); }

static inline VEC nll_gbls_vec(const VEC R2, float b) {
    const float a = 0.5f;
    const float ln_beta = fast_lbetaf_simd(a, b);
    const float split_threshold = (b + 1.0f) / (a + b + 2.0f);

    VEC out;
    for (int lane = 0; lane < VEC_LEN; ++lane) {
        float r2 = R2.values[lane];
        if (!float_is_finite_bits(r2) || r2 <= 0.0f) {
            out.values[lane] = 0.0f;
            continue;
        }
        if (r2 > AOV_R2_MAX) r2 = AOV_R2_MAX;

        float y = 1.0f - r2;
        float x = r2;

        if (y < split_threshold) {
            /* Tail evaluation: Iy(b, a) */
            float ln_front = b * logf(y) + a * logf(x) - logf(b) - ln_beta;
            float f = 1.0f, c = 1.0f, d_cf = 1.0f - (a + b) * y / (b + 1.0f);
            if (fabsf(d_cf) < 1.0e-30f) d_cf = 1.0e-30f;
            d_cf = 1.0f / d_cf;
            f *= d_cf;

            for (int m = 1; m <= 15; ++m) {
                float dm = (float)m;
                float num_even = dm * (a - dm) * y / ((b + 2.0f * dm - 1.0f) * (b + 2.0f * dm));
                d_cf = 1.0f + num_even * d_cf;
                if (fabsf(d_cf) < 1.0e-30f) d_cf = 1.0e-30f;
                c = 1.0f + num_even / c;
                if (fabsf(c) < 1.0e-30f) c = 1.0e-30f;
                d_cf = 1.0f / d_cf;
                f *= (c * d_cf);

                float num_odd = -(b + dm) * (a + b + dm) * y / ((b + 2.0f * dm) * (b + 2.0f * dm + 1.0f));
                d_cf = 1.0f + num_odd * d_cf;
                if (fabsf(d_cf) < 1.0e-30f) d_cf = 1.0e-30f;
                c = 1.0f + num_odd / c;
                if (fabsf(c) < 1.0e-30f) c = 1.0e-30f;
                d_cf = 1.0f / d_cf;
                float delta = c * d_cf;
                f *= delta;
                if (fabsf(delta - 1.0f) < 1.0e-7f) break;
            }
            float val = -(ln_front + logf(f));
            out.values[lane] = (val > 0.0f && float_is_finite_bits(val)) ? val : 0.0f;
        } else {
            /* Small R2 evaluation: Ix(a, b) */
            float front = expf(a * logf(x) + b * logf(y) - ln_beta) / a;
            float f = 1.0f, c = 1.0f, d_cf = 1.0f - (a + b) * x / (a + 1.0f);
            if (fabsf(d_cf) < 1.0e-30f) d_cf = 1.0e-30f;
            d_cf = 1.0f / d_cf;
            f *= d_cf;

            for (int m = 1; m <= 15; ++m) {
                float dm = (float)m;
                float num_even = dm * (b - dm) * x / ((a + 2.0f * dm - 1.0f) * (a + 2.0f * dm));
                d_cf = 1.0f + num_even * d_cf;
                if (fabsf(d_cf) < 1.0e-30f) d_cf = 1.0e-30f;
                c = 1.0f + num_even / c;
                if (fabsf(c) < 1.0e-30f) c = 1.0e-30f;
                d_cf = 1.0f / d_cf;
                f *= (c * d_cf);

                float num_odd = -(a + dm) * (a + b + dm) * x / ((a + 2.0f * dm) * (a + 2.0f * dm + 1.0f));
                d_cf = 1.0f + num_odd * d_cf;
                if (fabsf(d_cf) < 1.0e-30f) d_cf = 1.0e-30f;
                c = 1.0f + num_odd / c;
                if (fabsf(c) < 1.0e-30f) c = 1.0e-30f;
                d_cf = 1.0f / d_cf;
                float delta = c * d_cf;
                f *= delta;
                if (fabsf(delta - 1.0f) < 1.0e-7f) break;
            }
            float Ix = front * f;
            float val = (Ix <= 0.0f) ? 0.0f : -log1pf(-Ix);
            out.values[lane] = (val > 0.0f && float_is_finite_bits(val)) ? val : 0.0f;
        }
    }
    return out;
}

static inline DVEC nll_gbls_dvec(const DVEC R2, double b) {
    const double a = 0.5;
    const double ln_beta = fast_lbeta_simd(a, b);
    const double split_threshold = (b + 1.0) / (a + b + 2.0);

    DVEC out;
    for (int lane = 0; lane < DVEC_LEN; ++lane) {
        double r2 = R2.values[lane];
        if (!double_is_finite_bits(r2) || r2 <= 0.0) {
            out.values[lane] = 0.0;
            continue;
        }
        if (r2 > (double)AOV_R2_MAX) r2 = (double)AOV_R2_MAX;

        double y = 1.0 - r2;
        double x = r2;

        if (y < split_threshold) {
            /* Tail evaluation: Iy(b, a) */
            double ln_front = b * log(y) + a * log(x) - log(b) - ln_beta;
            double f = 1.0, c = 1.0, d_cf = 1.0 - (a + b) * y / (b + 1.0);
            if (fabs(d_cf) < 1.0e-30) d_cf = 1.0e-30;
            d_cf = 1.0 / d_cf;
            f *= d_cf;

            for (int m = 1; m <= 25; ++m) {
                double dm = (double)m;
                double num_even = dm * (a - dm) * y / ((b + 2.0 * dm - 1.0) * (b + 2.0 * dm));
                d_cf = 1.0 + num_even * d_cf;
                if (fabs(d_cf) < 1.0e-30) d_cf = 1.0e-30;
                c = 1.0 + num_even / c;
                if (fabs(c) < 1.0e-30) c = 1.0e-30;
                d_cf = 1.0 / d_cf;
                f *= (c * d_cf);

                double num_odd = -(b + dm) * (a + b + dm) * y / ((b + 2.0 * dm) * (b + 2.0 * dm + 1.0));
                d_cf = 1.0 + num_odd * d_cf;
                if (fabs(d_cf) < 1.0e-30) d_cf = 1.0e-30;
                c = 1.0 + num_odd / c;
                if (fabs(c) < 1.0e-30) c = 1.0e-30;
                d_cf = 1.0 / d_cf;
                float delta = c * d_cf;
                f *= delta;
                if (fabs(delta - 1.0) < 1.0e-15) break;
            }
            double val = -(ln_front + log(f));
            out.values[lane] = (val > 0.0 && double_is_finite_bits(val)) ? val : 0.0;
        } else {
            /* Small R2 evaluation: Ix(a, b) */
            double front = exp(a * log(x) + b * log(y) - ln_beta) / a;
            double f = 1.0, c = 1.0, d_cf = 1.0 - (a + b) * x / (a + 1.0);
            if (fabs(d_cf) < 1.0e-30) d_cf = 1.0e-30;
            d_cf = 1.0 / d_cf;
            f *= d_cf;

            for (int m = 1; m <= 25; ++m) {
                double dm = (double)m;
                double num_even = dm * (b - dm) * x / ((a + 2.0 * dm - 1.0) * (a + 2.0 * dm));
                d_cf = 1.0 + num_even * d_cf;
                if (fabs(d_cf) < 1.0e-30) d_cf = 1.0e-30;
                c = 1.0 + num_even / c;
                if (fabs(c) < 1.0e-30) c = 1.0e-30;
                d_cf = 1.0 / d_cf;
                f *= (c * d_cf);

                double num_odd = -(a + dm) * (a + b + dm) * x / ((a + 2.0 * dm) * (a + 2.0 * dm + 1.0));
                d_cf = 1.0 + num_odd * d_cf;
                if (fabs(d_cf) < 1.0e-30) d_cf = 1.0e-30;
                c = 1.0 + num_odd / c;
                if (fabs(c) < 1.0e-30) c = 1.0e-30;
                d_cf = 1.0 / d_cf;
                double delta = c * d_cf;
                f *= delta;
                if (fabs(delta - 1.0) < 1.0e-15) break;
            }
            double Ix = front * f;
            double val = (Ix <= 0.0) ? 0.0 : -log1p(-Ix);
            out.values[lane] = (val > 0.0 && double_is_finite_bits(val)) ? val : 0.0;
        }
    }
    return out;
}

/* ========================================================================= */
/* Unified lnFAP & Spectrum Batch Converter via Pure Vectorized Routines    */
/* ========================================================================= */

static inline double lnFAP_fast(int dK, int dH, double R2, int N) {
    if (!double_is_finite_bits(R2) || R2 <= 0.0) return 0.0;
    if (R2 >= 1.0) return INFINITY;

    int d = dK - dH;
    int Nk = N - dK;
    if (d < 1 || Nk < 1) return 0.0;

    if ((d % 2) == 0) {
        int degree = d / 2;
        if (degree <= 8) {
            VEC in = SET_VEC((float)R2);
            VEC out = nll_exact_pochhammer_vec(degree, in, 0.5f * (float)Nk);
            return (double)out.values[0];
        } else {
            DVEC in = SET_DVEC(R2);
            DVEC out = nll_exact_pochhammer_dvec(degree, in, 0.5 * (double)Nk);
            return out.values[0];
        }
    }

    /* GBLS / odd d (e.g. d = 1) */
    VEC in = SET_VEC((float)R2);
    VEC out = nll_gbls_vec(in, 0.5f * (float)Nk);
    return (double)out.values[0];
}

static inline float lnFAP_fast_f32(int dK, int dH, float R2, int N) {
    if (!float_is_finite_bits(R2) || R2 <= 0.0f) return 0.0f;
    if (R2 >= 1.0f) return INFINITY;

    int d = dK - dH;
    int Nk = N - dK;
    if (d < 1 || Nk < 1) return 0.0f;

    if ((d % 2) == 0) {
        int degree = d / 2;
        if (degree <= 8) {
            VEC in = SET_VEC(R2);
            VEC out = nll_exact_pochhammer_vec(degree, in, 0.5f * (float)Nk);
            return out.values[0];
        } else {
            DVEC in = SET_DVEC((double)R2);
            DVEC out = nll_exact_pochhammer_dvec(degree, in, 0.5 * (double)Nk);
            return (float)out.values[0];
        }
    }

    VEC in = SET_VEC(R2);
    VEC out = nll_gbls_vec(in, 0.5f * (float)Nk);
    return out.values[0];
}

/* Fast batch converter for full spectrum arrays (used by AoV and GBLS) */
static inline void nll_convert_spectrum_batch(const float *r2_in, float *nll_out, size_t count, int degree, int n_eff) {
    if (count == 0 || n_eff <= 2 * degree + 1) return;
    float b_f = 0.5f * (float)(n_eff - (2 * degree + 1));
    double b_d = 0.5 * (double)(n_eff - (2 * degree + 1));
    size_t i = 0;

    if (degree <= 8) {
        for (; i + VEC_LEN <= count; i += VEC_LEN) {
            VEC r2_chunk;
            memcpy(r2_chunk.values, r2_in + i, sizeof(float) * VEC_LEN);
            VEC out = nll_exact_pochhammer_vec(degree, r2_chunk, b_f);
            memcpy(nll_out + i, out.values, sizeof(float) * VEC_LEN);
        }
    } else {
        for (; i + DVEC_LEN <= count; i += DVEC_LEN) {
            DVEC r2_chunk;
            for (int lane = 0; lane < DVEC_LEN; ++lane) r2_chunk.values[lane] = (double)r2_in[i + lane];
            DVEC out = nll_exact_pochhammer_dvec(degree, r2_chunk, b_d);
            for (int lane = 0; lane < DVEC_LEN; ++lane) nll_out[i + lane] = (float)out.values[lane];
        }
    }

    if (i < count) {
        size_t rem = count - i;
        if (degree <= 8) {
            VEC r2_pad = {0};
            for (size_t k = 0; k < rem; ++k) r2_pad.values[k] = r2_in[i + k];
            VEC out = nll_exact_pochhammer_vec(degree, r2_pad, b_f);
            for (size_t k = 0; k < rem; ++k) nll_out[i + k] = out.values[k];
        } else {
            DVEC d_pad = {0};
            for (size_t k = 0; k < rem; ++k) d_pad.values[k] = (double)r2_in[i + k];
            DVEC out = nll_exact_pochhammer_dvec(degree, d_pad, b_d);
            for (size_t k = 0; k < rem; ++k) nll_out[i + k] = (float)out.values[k];
        }
    }
}

#ifndef lnFAP
#    define lnFAP(dK, dH, R2, N) lnFAP_fast((dK), (dH), (double)(R2), (N))
#endif

#endif  // SIMD_H
