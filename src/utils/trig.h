#ifndef TRIG_H
#define TRIG_H

#include <stdint.h>

#include "simd.h"

/* =========================================================================
 * Scalar quadrant-reduced sin/cos kernels.
 *
 * The phase argument is a float number of cycles; the fractional part is
 * taken with an (int) cast, so callers must keep |x| < 2^31.  Reference
 * consumers are expected to reduce large phases in double precision first
 * (see aov_fold_phase_dvec() and the twiddle ladder) and only then convert
 * to float.
 * ========================================================================= */

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

/* =========================================================================
 * Vector kernels and the double->float phase fold used by the twiddle ladder.
 *
 * The ladder resynchronises by folding an exact double phase into its
 * fractional part and regenerating cos/sin with the same quadrant-reduced
 * polynomials as above, but VEC-wide: two DVEC lanes are converted to two
 * half-wide float vectors and concatenated into one VEC.
 *
 * The fold truncates towards zero (the fractional part keeps the sign and the
 * vector reducers below apply the scalar "if (f < 0) f += 1" rule).
 *
 * By default the truncation uses a DVEC-sized 32-bit integer vector: SSE2 and
 * AVX2 have a vector double->int32 conversion (cvttpd2dq) but only a scalar
 * double->int64 one, so the 32-bit form is shorter and branch-free (measured
 * on the fold loop: 13 vs 19 instructions on SSE2, 15 vs 27 on AVX2, neutral on
 * AVX-512).  The 32-bit truncation assumes |phase| < 2^31, which is exactly the
 * contract of the scalar kernels above.  Build with -DLADDER_FOLD_INT64=1 for
 * wider phases; that variant is bit-identical for canonical phases.
 * ========================================================================= */

/* Everything below relies on GNU vector extensions plus __builtin_convertvector,
 * exactly like simd.h already does for the rest of the SIMD layer. */

#ifndef LADDER_FOLD_INT64
#    define LADDER_FOLD_INT64 0
#endif

/* One lane per double in a DVEC. */
#if LADDER_FOLD_INT64
typedef int64_t trig_fold_int __attribute__((vector_size(VEC_BYTES)));
#else
typedef int32_t trig_fold_int __attribute__((vector_size(VEC_BYTES / 2)));
#endif
/* Half a VEC: exactly DVEC_LEN floats. */
typedef float trig_half_vecf __attribute__((vector_size(VEC_BYTES / 2)));

typedef union {
    vecf_data data;
    trig_half_vecf halves[2];
} trig_vecf_pair;

/* Exact fractional part of a double phase: phase - trunc(phase), |result| < 1. */
static inline DVEC aov_fold_phase_dvec(DVEC phase) {
    DVEC out;
    trig_fold_int trunc_int = __builtin_convertvector(phase.data, trig_fold_int);
    out.data = phase.data - __builtin_convertvector(trunc_int, vecd_data);
    return out;
}

/* Concatenate two DVEC fractional phases into one VEC (lane order preserved). */
static inline VEC aov_dvec_pair_to_vecf(DVEC lo, DVEC hi) {
    trig_vecf_pair pair;
    pair.halves[0] = __builtin_convertvector(lo.data, trig_half_vecf);
    pair.halves[1] = __builtin_convertvector(hi.data, trig_half_vecf);
    return (VEC){.data = pair.data};
}

/* Vector counterparts of sin2pif_tls()/cos2pif_tls() (same constants and
 * reduction).  Comparison results are -1/0 masks, which doubles as the sign. */
static inline VEC sin2pif_vec(VEC x) {
    VEC f = {.data = x.data - __builtin_convertvector(__builtin_convertvector(x.data, veci_data), vecf_data)};
    veci_data neg = f.data < (vecf_data){};
    f.data = f.data - __builtin_convertvector(neg, vecf_data);
    veci_data big = f.data >= ((vecf_data){} + 0.5f);
    veci_data sign = ((veci_data){} + 1) + 2 * big;
    f.data = f.data + __builtin_convertvector(big, vecf_data) * 0.5f;
    veci_data fold = f.data > ((vecf_data){} + 0.25f);
    f.data = f.data + __builtin_convertvector(fold, vecf_data) * (2.0f * f.data - 0.5f);
    vecf_data f2 = f.data * f.data;
    vecf_data p = f2 * 39.536706065730207835108712734262f - 76.549782293595742666226937116116f;
    p = p * f2 + 81.601004073261773523492199897936f;
    p = p * f2 - 41.341655031416278077153126232486f;
    p = p * f2 + 6.2831851600894774430188071795666f;
    p = p * f.data;
    return (VEC){.data = p * __builtin_convertvector(sign, vecf_data)};
}

static inline VEC cos2pif_vec(VEC x) {
    VEC f = {.data = x.data - __builtin_convertvector(__builtin_convertvector(x.data, veci_data), vecf_data)};
    veci_data neg = f.data < (vecf_data){};
    f.data = f.data - __builtin_convertvector(neg, vecf_data);
    veci_data over = f.data > ((vecf_data){} + 0.5f);
    f.data = f.data + __builtin_convertvector(over, vecf_data) * (2.0f * f.data - 1.0f);
    veci_data fold = f.data > ((vecf_data){} + 0.25f);
    veci_data sign = ((veci_data){} + 1) + 2 * fold;
    f.data = f.data + __builtin_convertvector(fold, vecf_data) * (2.0f * f.data - 0.5f);
    vecf_data f2 = f.data * f.data;
    vecf_data p = f2 * 56.242380464873243259663276802701f - 85.240330322699427859509454517828f;
    p = p * f2 + 64.934590626780991246193352727536f;
    p = p * f2 - 19.739171434702393618770795066531f;
    p = p * f2 + 0.99999995346667013630639784578184f;
    return (VEC){.data = p * __builtin_convertvector(sign, vecf_data)};
}

#endif  // TRIG_H
