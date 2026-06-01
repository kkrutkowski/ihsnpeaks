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
