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
        __m256 data;
        uint32_t values[8];
    } m256i_union;
    #define VEC m256_union
    #define IVEC m256i_union
#endif


#endif // SIMD_H
