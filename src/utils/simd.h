#ifndef SIMD_H
#define SIMD_H
#include <stdint.h>
#include <math.h>

// Type definitions
#ifdef __AVX512F__
        #define SET_VEC(val) ((m512_union){.values={val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val}})
        #define SET_IVEC(val) ((m512i_union){.values={val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val}})
        #define SET_DVEC(val) ((m512_union){.values={val, val, val, val, val, val, val, val}})

    #include <immintrin.h>
    typedef union {__m512 data; float values[16];} m512_union;
    typedef union {__m512d data; double values[8];} m512d_union;
    typedef union {__m512i data; int32_t values[16];} m512i_union;

    typedef m512_union VEC;
    typedef m512d_union DVEC;
    typedef m512i_union IVEC;
#else
    #ifdef __AVX__
        #define SET_VEC(val) ((m256_union){.values={val, val, val, val, val, val, val, val}})
        #define SET_IVEC(val) ((m256i_union){.values={val, val, val, val, val, val, val, val}})
        #define SET_DVEC(val) ((m256_union){.values={val, val, val, val}})

        #include <immintrin.h>
        typedef union {__m256 data; float values[8];} m256_union;
        typedef union {__m256d data; double values[4];} m256d_union;
        typedef union {__m256i data; int32_t values[8];} m256i_union;

        typedef m256_union VEC;
        typedef m256d_union DVEC;
        typedef m256i_union IVEC;
    #else
        // GNU Vector extensions for 256-bit vectors
        typedef float v8sf __attribute__ ((vector_size (32)));
        typedef double v4df __attribute__ ((vector_size (32)));
        typedef int32_t v8si __attribute__ ((vector_size (32)));

        #define SET_VEC(val) ((v8sf){.values={val, val, val, val, val, val, val, val}})
        #define SET_IVEC(val) ((v8si){.values={val, val, val, val, val, val, val, val}})
        #define SET_DVEC(val) ((v4df){.values={val, val, val, val}})

        typedef union {v8sf data; float values[8];} m256_union;
        typedef union {v4df data; double values[4];} m256d_union;
        typedef union {v8si data; int32_t values[8];} m256i_union;

        typedef m256_union VEC;
        typedef m256d_union DVEC;
        typedef m256i_union IVEC;

        // Helper function for floor using __builtin_convertvector
        static inline v8sf vec_floor_ps(v8sf x) {
            v8si truncated = __builtin_convertvector(x, v8si);
            return __builtin_convertvector(truncated, v8sf);
        }
    #endif
#endif

#define VEC_LEN (sizeof(VEC) / sizeof(float))
#define DVEC_LEN (sizeof(DVEC) / sizeof(double))
#define IVEC_LEN (sizeof(IVEC) / sizeof(int32_t))

// Definitions of functions
#ifdef __AVX512F__
    static inline VEC sin_2pi_poly_ps(const VEC x) {
        constexpr VEC c[4] = {SET_VEC(6.2831676e+0f), SET_VEC(-4.1337518e+1f), SET_VEC(8.1351678e+1f), SET_VEC(-7.1087358e+1f)};
        const __m512 x2 = _mm512_mul_ps(x.data, x.data);
        VEC result;
        result.data = _mm512_fmadd_ps(c[3].data, x2, c[2].data);
        result.data = _mm512_fmadd_ps(result.data, x2, c[1].data);
        result.data = _mm512_fmadd_ps(result.data, x2, c[0].data);
        result.data = _mm512_mul_ps(result.data, x.data);
        return result;
    }

    static inline VEC sin_2pi_ps(const VEC angle) {
        constexpr VEC c[4] = {SET_VEC(0.25f), SET_VEC(0.5f), SET_VEC(0.75f), SET_VEC(1.0f)};
        const __m512 AVX512_SIGNMASK_PS = _mm512_castsi512_ps(_mm512_set1_epi32(0x80000000));
        VEC sinangle;
        sinangle.data = _mm512_sub_ps(angle.data, _mm512_floor_ps(angle.data));
        const __m512 angle_orig = sinangle.data;

        __mmask16 mask0 = _mm512_cmp_ps_mask(angle_orig, c[0].data, _CMP_GE_OQ);
        sinangle.data = _mm512_mask_sub_ps(sinangle.data, mask0, c[1].data, angle_orig);

        __mmask16 mask1 = _mm512_cmp_ps_mask(angle_orig, c[1].data, _CMP_GE_OQ);
        sinangle.data = _mm512_mask_xor_ps(sinangle.data, mask1, sinangle.data, AVX512_SIGNMASK_PS);

        __mmask16 mask2 = _mm512_cmp_ps_mask(angle_orig, c[2].data, _CMP_GE_OQ);
        sinangle.data = _mm512_mask_sub_ps(sinangle.data, mask2, c[3].data, angle_orig);

        VEC result = sin_2pi_poly_ps(sinangle);
        result.data = _mm512_mask_xor_ps(result.data, mask1, result.data, AVX512_SIGNMASK_PS);
        return result;
    }

    static inline VEC generateWeights(const float dst) {
        constexpr VEC c[4] = {SET_VEC(0.16666666f), SET_VEC(0.5f), SET_VEC(3.0f), SET_VEC(M_PI * M_PI)};
        constexpr VEC DST = (m512_union){{7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, -1.0f, -2.0f, -3.0f, -4.0f, -5.0f, -6.0f, -7.0f, -8.0f}};

        VEC denom = { .data = _mm512_add_ps(DST.data, _mm512_set1_ps(dst)) };
        VEC temp1 = { .data = _mm512_mul_ps(c[0].data, denom.data) };
        VEC temp2 = { .data = _mm512_mul_ps(c[1].data, denom.data) };
        VEC num = { .data = _mm512_mul_ps(c[2].data, _mm512_mul_ps(sin_2pi_ps(temp1).data, sin_2pi_ps(temp2).data)) };
        denom.data = _mm512_mul_ps(_mm512_mul_ps(denom.data, denom.data), c[3].data);
        VEC weights = { .data = _mm512_div_ps(num.data, denom.data) };
        return weights;
    }
#else
    #ifdef __AVX__
        static inline VEC sin_2pi_poly_ps(const VEC x) {
            constexpr VEC c[4] = {SET_VEC(6.2831676e+0f), SET_VEC(-4.1337518e+1f), SET_VEC(8.1351678e+1f), SET_VEC(-7.1087358e+1f)};
            const __m256 x2 = _mm256_mul_ps(x.data, x.data);

            VEC result;
            #ifdef __FMA__
                result.data = _mm256_fmadd_ps(c[3].data, x2, c[2].data);
                result.data = _mm256_fmadd_ps(result.data, x2, c[1].data);
                result.data = _mm256_fmadd_ps(result.data, x2, c[0].data);
            #else
                result.data = _mm256_add_ps(_mm256_mul_ps(c[3].data, x2), c[2].data);
                result.data = _mm256_add_ps(_mm256_mul_ps(result.data, x2), c[1].data);
                result.data = _mm256_add_ps(_mm256_mul_ps(result.data, x2), c[0].data);
            #endif
            result.data = _mm256_mul_ps(result.data, x.data);
            return result;
        }

        static inline VEC sin_2pi_ps(const VEC angle) {
            constexpr VEC c[4] = {SET_VEC(0.25f), SET_VEC(0.5f), SET_VEC(0.75f), SET_VEC(1.0f)};
            const __m256 AVX_SIGNMASK_PS = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));

            VEC sinangle;
            sinangle.data = _mm256_sub_ps(angle.data, _mm256_floor_ps(angle.data));
            const __m256 angle_orig = sinangle.data;

            sinangle.data = _mm256_xor_ps(sinangle.data, _mm256_and_ps(_mm256_cmp_ps(angle_orig, c[0].data, _CMP_GE_OQ), _mm256_xor_ps(sinangle.data, _mm256_sub_ps(c[1].data, angle_orig))));
            sinangle.data = _mm256_xor_ps(sinangle.data, _mm256_and_ps(_mm256_cmp_ps(angle_orig, c[1].data, _CMP_GE_OQ), AVX_SIGNMASK_PS));
            sinangle.data = _mm256_xor_ps(sinangle.data, _mm256_and_ps(_mm256_cmp_ps(angle_orig, c[2].data, _CMP_GE_OQ), _mm256_xor_ps(sinangle.data, _mm256_sub_ps(c[3].data, angle_orig))));

            VEC result = sin_2pi_poly_ps(sinangle);
            result.data = _mm256_xor_ps(result.data, _mm256_and_ps(_mm256_cmp_ps(angle_orig, c[1].data, _CMP_GE_OQ), AVX_SIGNMASK_PS));
            return result;
        }

        static inline void generateWeights(const float dst, VEC *h1, VEC *h2) {
            constexpr VEC c[4] = {SET_VEC(0.16666666f), SET_VEC(0.5f), SET_VEC(3.0f), SET_VEC(M_PI * M_PI)};
            constexpr VEC DST[2] = {(m256_union){{7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f}}, (m256_union){{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f}}};

            VEC denom[2];
            denom[0].data = _mm256_add_ps(DST[0].data, _mm256_set1_ps(dst));
            denom[1].data = _mm256_sub_ps(DST[1].data, _mm256_set1_ps(dst));

            VEC num[2];
            num[0].data = _mm256_mul_ps(c[2].data, _mm256_mul_ps(sin_2pi_ps((VEC){.data = _mm256_mul_ps(c[0].data, denom[0].data)}).data, sin_2pi_ps((VEC){.data = _mm256_mul_ps(c[1].data, denom[0].data)}).data));
            num[1].data = _mm256_mul_ps(c[2].data, _mm256_mul_ps(sin_2pi_ps((VEC){.data = _mm256_mul_ps(c[0].data, denom[1].data)}).data, sin_2pi_ps((VEC){.data = _mm256_mul_ps(c[1].data, denom[1].data)}).data));

            denom[0].data = _mm256_mul_ps(_mm256_mul_ps(denom[0].data, denom[0].data), c[3].data);
            denom[1].data = _mm256_mul_ps(_mm256_mul_ps(denom[1].data, denom[1].data), c[3].data);

            h1->data = _mm256_div_ps(num[0].data, denom[0].data);
            h2->data = _mm256_div_ps(num[1].data, denom[1].data);
        }
    #else
        // GNU Vector extensions implementation using built-in operators
        static inline VEC sin_2pi_poly_ps(const VEC x) {
            constexpr VEC c[4] = {SET_VEC(6.2831676e+0f), SET_VEC(-4.1337518e+1f), SET_VEC(8.1351678e+1f), SET_VEC(-7.1087358e+1f)};
            const v8sf x2 = x.data * x.data;

            VEC result;
            result.data = c[3].data * x2 + c[2].data;
            result.data = result.data * x2 + c[1].data;
            result.data = result.data * x2 + c[0].data;
            result.data = result.data * x.data;
            return result;
        }

        static inline VEC sin_2pi_ps(const VEC angle) {
            constexpr VEC c[4] = {SET_VEC(0.25f), SET_VEC(0.5f), SET_VEC(0.75f), SET_VEC(1.0f)};
            const v8sf AVX_SIGNMASK_PS = (v8sf){-0.0f, -0.0f, -0.0f, -0.0f, -0.0f, -0.0f, -0.0f, -0.0f};

            VEC sinangle;
            sinangle.data = angle.data - vec_floor_ps(angle.data);
            const v8sf angle_orig = sinangle.data;

            // Use built-in comparison operators that return mask vectors
            v8sf mask0 = angle_orig >= c[0].data;
            v8sf mask1 = angle_orig >= c[1].data;
            v8sf mask2 = angle_orig >= c[2].data;

            sinangle.data = sinangle.data ^ (mask0 & (sinangle.data ^ (c[1].data - angle_orig)));
            sinangle.data = sinangle.data ^ (mask1 & AVX_SIGNMASK_PS);
            sinangle.data = sinangle.data ^ (mask2 & (sinangle.data ^ (c[3].data - angle_orig)));

            VEC result = sin_2pi_poly_ps(sinangle);
            result.data = result.data ^ (mask1 & AVX_SIGNMASK_PS);
            return result;
        }

        static inline void generateWeights(const float dst, VEC *h1, VEC *h2) {
            constexpr VEC c[4] = {SET_VEC(0.16666666f), SET_VEC(0.5f), SET_VEC(3.0f), SET_VEC(M_PI * M_PI)};
            constexpr VEC DST[2] = {(m256_union){{7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f}}, (m256_union){{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f}}};

            VEC denom[2];
            const v8sf dst_vec = (v8sf){dst, dst, dst, dst, dst, dst, dst, dst};
            denom[0].data = DST[0].data + dst_vec;
            denom[1].data = DST[1].data - dst_vec;

            VEC num[2];
            VEC temp0 = { .data = c[0].data * denom[0].data };
            VEC temp1 = { .data = c[1].data * denom[0].data };
            num[0].data = c[2].data * (sin_2pi_ps(temp0).data * sin_2pi_ps(temp1).data);

            temp0.data = c[0].data * denom[1].data;
            temp1.data = c[1].data * denom[1].data;
            num[1].data = c[2].data * (sin_2pi_ps(temp0).data * sin_2pi_ps(temp1).data);

            denom[0].data = (denom[0].data * denom[0].data) * c[3].data;
            denom[1].data = (denom[1].data * denom[1].data) * c[3].data;

            h1->data = num[0].data / denom[0].data;
            h2->data = num[1].data / denom[1].data;
        }
    #endif
#endif

// Natural logarithm implementation using GNU vector extensions
static inline VEC ln_ps(const VEC x) {
    // Constants for frexp
    constexpr IVEC exp_mask = SET_IVEC(0x7F800000);  // Mask for exponent bits
    constexpr IVEC mant_mask = SET_IVEC(0x807FFFFF); // Mask for mantissa (clear exp)
    constexpr VEC one_vec = SET_VEC(1.0f);           // Value 1.0
    constexpr VEC ln2_vec = SET_VEC(0.69314718f);    // ln(2)
    constexpr IVEC exp_bias = SET_IVEC(127);         // Exponent bias
    constexpr IVEC set_1 = SET_IVEC(0x3F800000);     // Bit pattern for 1.0

    // Coefficients for ln(1 + x) on [0, 1)
    constexpr VEC c[6] = {
        SET_VEC(-1.936759742e0f), SET_VEC(3.514087297e0f),
        SET_VEC(-2.440029763e0f), SET_VEC(1.116090027e0f),
        SET_VEC(-2.83826848e-1f), SET_VEC(3.04490045e-2f)
    };

    IVEC x_bits, exp_bits, unbiased_exp, mant_bits; VEC mant_vec, exp_float, exp_ln2, x_minus_one, ln_result; // Initialize the temporary vectors

    // Convert float vector to int vector for bit manipulation
    x_bits.data = (typeof(x_bits.data))x.data;

    // frexp: Extract exponent and mantissa
    exp_bits.data = (x_bits.data >> 23) & SET_IVEC(0xFF).data; // Extract exponent bits
    unbiased_exp.data = exp_bits.data - exp_bias.data;       // Unbias exponent

    // Extract mantissa: clear exponent bits and set them to 0x3F800000 (1.0)
    mant_bits.data = (x_bits.data & mant_mask.data) | set_1.data;
    mant_vec.data = (typeof(mant_vec.data))mant_bits.data;

    // Compute e * ln(2) - convert integer exponent to float
    #ifdef __AVX512F__
        exp_float.data = _mm512_cvtepi32_ps(unbiased_exp.data);
    #elif defined(__AVX__)
        exp_float.data = _mm256_cvtepi32_ps(unbiased_exp.data);
    #else
        exp_float.data = __builtin_convertvector(unbiased_exp.data, typeof(exp_float.data));
   #endif
    exp_ln2.data = exp_float.data * ln2_vec.data;

    // Horner's method for polynomial evaluation
    VEC ln_mant = c[5];
    ln_mant.data = ln_mant.data * mant_vec.data + c[4].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[3].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[2].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[1].data;
    ln_mant.data = ln_mant.data * mant_vec.data + c[0].data;

    // Combine results: ln(x) = e * ln(2) + ln(mantissa)
    ln_result.data = exp_ln2.data + ln_mant.data;

    return ln_result;
}

#endif // SIMD_H
