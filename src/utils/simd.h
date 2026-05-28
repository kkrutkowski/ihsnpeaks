#ifndef SIMD_H
#define SIMD_H

#include <math.h>
#include <stdint.h>

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

typedef float vecf_data __attribute__((vector_size(VEC_BYTES)));
typedef double vecd_data __attribute__((vector_size(VEC_BYTES)));
typedef int32_t veci_data __attribute__((vector_size(VEC_BYTES)));

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

#define SET_VEC(val) ((VEC){.data = (vecf_data){} + (float)(val)})
#define SET_DVEC(val) ((DVEC){.data = (vecd_data){} + (double)(val)})
#define SET_IVEC(val) ((IVEC){.data = (veci_data){} + (int32_t)(val)})

static inline VEC vec_blend(const veci_data mask, const VEC when_true, const VEC when_false) {
    union {
        vecf_data f;
        veci_data i;
    } t = {.f = when_true.data}, f = {.f = when_false.data}, out;
    out.i = f.i ^ (mask & (f.i ^ t.i));
    return (VEC){.data = out.f};
}

static inline VEC vec_trunc_ps(const VEC x) {
    IVEC truncated;
    truncated.data = __builtin_convertvector(x.data, veci_data);
    return (VEC){.data = __builtin_convertvector(truncated.data, vecf_data)};
}

static inline float sin2pif_tls(float x) {
    float f = x - (float)((int)x);
    if (f < 0.0f) f += 1.0f;
    float sign = 1.0f;
    if (f >= 0.5f) {
        sign = -1.0f;
        f -= 0.5f;
    }
    if (f > 0.25f) f = 0.5f - f;
    float f2 = f * f;
    float p = f2 * 39.536706065730207835108712734262f - 76.549782293595742666226937116116f;
    p = p * f2 + 81.601004073261773523492199897936f;
    p = p * f2 - 41.341655031416278077153126232486f;
    p = p * f2 + 6.2831851600894774430188071795666f;
    p *= f;
    return p * sign;
}

static inline float cos2pif_tls(float x) {
    float f = x - (float)((int)x);
    if (f < 0.0f) f += 1.0f;
    if (f > 0.5f) f = 1.0f - f;
    float sign = 1.0f;
    if (f > 0.25f) {
        sign = -1.0f;
        f = 0.5f - f;
    }
    float f2 = f * f;
    float p = f2 * 56.242380464873243259663276802701f - 85.240330322699427859509454517828f;
    p = p * f2 + 64.934590626780991246193352727536f;
    p = p * f2 - 19.739171434702393618770795066531f;
    p = p * f2 + 0.99999995346667013630639784578184f;
    return p * sign;
}

static inline VEC sin_2pi_ps(const VEC x) {
    const VEC zero = SET_VEC(0.0f);
    const VEC one = SET_VEC(1.0f);
    const VEC half = SET_VEC(0.5f);
    const VEC quarter = SET_VEC(0.25f);
    const VEC minus_one = SET_VEC(-1.0f);

    VEC f = {.data = x.data - vec_trunc_ps(x).data};
    f = vec_blend(f.data < zero.data, (VEC){.data = f.data + one.data}, f);

    VEC sign = one;
    veci_data ge_half = f.data >= half.data;
    sign = vec_blend(ge_half, minus_one, sign);
    f = vec_blend(ge_half, (VEC){.data = f.data - half.data}, f);
    f = vec_blend(f.data > quarter.data, (VEC){.data = half.data - f.data}, f);

    VEC f2 = {.data = f.data * f.data};
    VEC p = {.data = f2.data * SET_VEC(39.536706065730207835108712734262f).data - SET_VEC(76.549782293595742666226937116116f).data};
    p.data = p.data * f2.data + SET_VEC(81.601004073261773523492199897936f).data;
    p.data = p.data * f2.data - SET_VEC(41.341655031416278077153126232486f).data;
    p.data = p.data * f2.data + SET_VEC(6.2831851600894774430188071795666f).data;
    return (VEC){.data = p.data * f.data * sign.data};
}

static inline VEC cos_2pi_ps(const VEC x) {
    const VEC zero = SET_VEC(0.0f);
    const VEC one = SET_VEC(1.0f);
    const VEC half = SET_VEC(0.5f);
    const VEC quarter = SET_VEC(0.25f);
    const VEC minus_one = SET_VEC(-1.0f);

    VEC f = {.data = x.data - vec_trunc_ps(x).data};
    f = vec_blend(f.data < zero.data, (VEC){.data = f.data + one.data}, f);
    f = vec_blend(f.data > half.data, (VEC){.data = one.data - f.data}, f);

    VEC sign = one;
    veci_data gt_quarter = f.data > quarter.data;
    sign = vec_blend(gt_quarter, minus_one, sign);
    f = vec_blend(gt_quarter, (VEC){.data = half.data - f.data}, f);

    VEC f2 = {.data = f.data * f.data};
    VEC p = {.data = f2.data * SET_VEC(56.242380464873243259663276802701f).data - SET_VEC(85.240330322699427859509454517828f).data};
    p.data = p.data * f2.data + SET_VEC(64.934590626780991246193352727536f).data;
    p.data = p.data * f2.data - SET_VEC(19.739171434702393618770795066531f).data;
    p.data = p.data * f2.data + SET_VEC(0.99999995346667013630639784578184f).data;
    return (VEC){.data = p.data * sign.data};
}

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

#endif  // SIMD_H
