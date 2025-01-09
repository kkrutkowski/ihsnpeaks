#ifndef SIMD_H
#define SIMD_H
#include <stdint.h>

#ifdef __AVX__
    #include <immintrin.h>
    typedef union {
        __m256 data;
        float values[8];
    } m256_union;

    typedef union {
        __m256i data;
        uint32_t values[8];
    } m256i_union;
    typedef m256_union VEC;
    typedef m256i_union IVEC;
#else
    typedef union {
        uint32_t data;
        float values[1];
    } scalar_union;

    typedef union {
        uint32_t data;
        uint32_t values[1];
    } scalari_union;
    typedef scalar_union VEC;
    typedef scalari_union IVEC;
#endif

#define VEC_LEN (sizeof(VEC) / sizeof(float))

#endif // SIMD_H
