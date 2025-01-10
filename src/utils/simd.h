#ifndef SIMD_H
#define SIMD_H
#include <stdint.h>

//type definitions
#ifdef __AVX__
    #define SET_VEC(val) ((m256_union){{val, val, val, val, val, val, val, val}})
    #define SET_IVEC(val) ((m256_unioni){{val, val, val, val, val, val, val, val}})
    #define SET_DVEC(val) ((m256_union){{val, val, val, val}})

    #include <immintrin.h>
    typedef union {__m256 data; float values[8];} m256_union;
    typedef union {__m256d data; double values[4];} m256d_union;
    typedef union {__m256i data;uint32_t values[8];} m256i_union;

    typedef m256_union VEC;
    typedef m256d_union DVEC;
    typedef m256i_union IVEC;
#else
    typedef union {uint32_t data; float values[1];} scalar_union;
    typedef union {uint64_t data; double values[1];} scalard_union;
    typedef union {uint32_t data; uint32_t values[1];} scalari_union;

    typedef scalar_union VEC;
    typedef scalard_union DVEC;
    typedef scalari_union IVEC;
#endif

#define VEC_LEN (sizeof(VEC) / sizeof(float))
#define DVEC_LEN (sizeof(DVEC) / sizeof(double))
#define IVEC_LEN (sizeof(IVEC) / sizeof(uint32_t))

//definitions of functions
#ifdef __AVX__
    static inline VEC sin_2pi_poly_ps(const VEC x) {
        constexpr VEC c[4] = {SET_VEC(6.2831676e+0f), SET_VEC(-4.1337518e+1f), SET_VEC(8.1351678e+1f), SET_VEC(-7.1087358e+1f)};
        const __m256 x2 = _mm256_mul_ps(x.data, x.data);

        // Using FMA instructions for odd powers
        VEC result;
        #ifdef __FMA__
            result.data = _mm256_fmadd_ps(c[3].data, x2, c[2].data);   // c3 + c4 * x^2
            result.data = _mm256_fmadd_ps(result.data, x2, c[1].data); // c2 + (c3 + c4 * x^2) * x^2
            result.data = _mm256_fmadd_ps(result.data, x2, c[0].data); // c1 + (c2 + (c3 + c4 * x^2) * x^2) * x^2
        #else
            result.data = _mm256_add_ps(_mm256_mul_ps(c[3].data, x2), c[2].data);
            result.data = _mm256_add_ps(_mm256_mul_ps(result.data, x2), c[1].data);
            result.data = _mm256_add_ps(_mm256_mul_ps(result.data, x2), c[0].data);
        #endif
        result.data = _mm256_mul_ps(result.data, x.data); // Multiply by x to apply odd powers
    return result;}

    static inline VEC sin_2pi_ps(const VEC angle) {
        constexpr VEC c[4] = {SET_VEC(0.25f), SET_VEC(0.5f), SET_VEC(0.75f), SET_VEC(1.0f)};
        const __m256 AVX_SIGNMASK_PS =  _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));

        VEC sinangle;
        sinangle.data = _mm256_sub_ps(angle.data, _mm256_floor_ps(angle.data));
        const __m256 angle_orig = sinangle.data;

        sinangle.data = _mm256_xor_ps(sinangle.data, _mm256_and_ps(_mm256_cmp_ps(angle_orig, c[0].data, _CMP_GE_OQ), _mm256_xor_ps(sinangle.data, _mm256_sub_ps(c[1].data, angle_orig))));
        sinangle.data = _mm256_xor_ps(sinangle.data,_mm256_and_ps(_mm256_cmp_ps(angle_orig, c[1].data, _CMP_GE_OQ), AVX_SIGNMASK_PS));
        sinangle.data = _mm256_xor_ps(sinangle.data, _mm256_and_ps(_mm256_cmp_ps(angle_orig, c[2].data, _CMP_GE_OQ),_mm256_xor_ps(sinangle.data, _mm256_sub_ps(c[3].data, angle_orig))));

        VEC result;
        result = sin_2pi_poly_ps(sinangle);
        result.data = _mm256_xor_ps(result.data, _mm256_and_ps(_mm256_cmp_ps(angle.data, c[1].data, _CMP_GE_OQ), AVX_SIGNMASK_PS));

    return result;
}

    static inline void generateWeights(const float dst, VEC *h1, VEC *h2){
        constexpr VEC c[4] = {SET_VEC(0.16666666f), SET_VEC(0.5f), SET_VEC(3.0f), SET_VEC(M_PI * M_PI)};
        constexpr VEC DST[2] = {(m256_union){7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f}, (m256_union){1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f}};

        VEC denom[2]; denom[0].data = _mm256_add_ps(DST[0].data, _mm256_set1_ps(dst)); denom[1].data = _mm256_sub_ps(DST[1].data, _mm256_set1_ps(dst));
        VEC num[2];
        num[0].data = c[2].data * sin_2pi_ps((m256_union)(c[0].data * denom[0].data)).data * sin_2pi_ps((m256_union)(c[1].data * denom[0].data)).data;
        num[1].data = c[2].data * sin_2pi_ps((m256_union)(c[0].data * denom[1].data)).data * sin_2pi_ps((m256_union)(c[1].data * denom[1].data)).data;

        denom[0].data = denom[0].data * denom[0].data * c[3].data; denom[1].data = denom[1].data * denom[1].data * c[3].data;

        h1->data = num[0].data / denom[0].data; h2->data = num[1].data / denom[1].data;
    return;}
#endif

#endif // SIMD_H
