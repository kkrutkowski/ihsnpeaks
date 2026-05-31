#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../utils/trig.h"
#include "nufft1.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#    define M_SQRT1_2 0.70710678118654752440
#endif

#ifndef MAX_TWIDDLE_REUSE
#    define MAX_TWIDDLE_REUSE 8
#endif
#if MAX_TWIDDLE_REUSE < 2 || (MAX_TWIDDLE_REUSE & (MAX_TWIDDLE_REUSE - 1)) != 0
#    error "MAX_TWIDDLE_REUSE must be a power of two greater than or equal to 2"
#endif

#define VECF_LEN (VEC_BYTES / 4)
typedef float VECF __attribute__((vector_size(VEC_BYTES)));
typedef uint32_t VECF_INT __attribute__((vector_size(VEC_BYTES)));
typedef int32_t VECF_SINT __attribute__((vector_size(VEC_BYTES)));

#define LOAD_VEC(ptr) (*(const VECF *)(ptr))
#define STORE_VEC(ptr, val) (*(VECF *)(ptr) = (val))

static inline void storeu_vecf(float *ptr, VECF val) { memcpy(ptr, &val, sizeof(val)); }

enum { PSWF_W = NUFFT1_PSWF_WIDTH, PSWF21_DEGREE = 8, PSWF43_DEGREE = 8, PSWF_BACKEND_43 = 1 };

static void *alloc_aligned_bytes(size_t size) {
    size_t alloc_size = size < 64 ? 64 : size;
    alloc_size = (alloc_size + 63) & ~(size_t)63;
    void *ptr = aligned_alloc(64, alloc_size);
    if (!ptr) return NULL;
    memset(ptr, 0, alloc_size);
    return ptr;
}

static float *alloc_aligned_float(size_t n) { return (float *)alloc_aligned_bytes(n * sizeof(float)); }

static inline size_t sz_max(size_t a, size_t b) { return a > b ? a : b; }

static inline int bitceil_int(unsigned int n) {
    if (n <= 1) return 1;
    unsigned int v = n - 1U;
    int shift = 0;
    while (v) {
        v >>= 1;
        ++shift;
    }
    return 1 << shift;
}

static int next_fast_len(int n) { return bitceil_int((unsigned int)n); }

static inline VECF vecf_exp_tls(VECF v) {
    VECF y = v * ((VECF){} + 1.4426950408889634f);
    VECF_SINT i = __builtin_convertvector(y, VECF_SINT);
    VECF i_f = __builtin_convertvector(i, VECF);
    VECF f = y - i_f;

    VECF_SINT is_neg = f < (VECF){};
    i = i + is_neg;
    f = f - __builtin_convertvector(is_neg, VECF);

    VECF p = f * 0.0f + 0.0018964611454333f;
    p = p * f + ((VECF){} + 0.0089428289841091f);
    p = p * f + ((VECF){} + 0.0558662463045207f);
    p = p * f + ((VECF){} + 0.2401397110907695f);
    p = p * f + ((VECF){} + 0.6931547524751674f);
    p = p * f + ((VECF){} + 0.9999998931108267f);

    union {
        VECF_SINT i;
        VECF f;
    } bitcast;
    bitcast.i = (i + ((VECF_SINT){} + 127)) << 23;
    return p * bitcast.f;
}

static inline float expf_remez_tls(float x) {
    VECF v = (VECF){};
    v[0] = x;
    return vecf_exp_tls(v)[0];
}

/* -------------------------------------------------------------------------
 * Single-precision nanoFFT
 * ------------------------------------------------------------------------- */

typedef struct nufft1_nanofft_plan {
    uint32_t n;
    float *twiddle_real;
    float *twiddle_imag;
} nufft1_nanofft_plan;

static void nufft1_nanofft_destroy_plan(nufft1_nanofft_plan *plan);

static inline bool is_power_of_two(uint32_t n) { return n > 0 && ((n & (n - 1U)) == 0); }
static inline uint32_t uint_min(uint32_t a, uint32_t b) { return a < b ? a : b; }
static inline uint32_t int_max(int32_t a, int32_t b) { return a > b ? (uint32_t)a : (uint32_t)b; }

#define LOG_BLOCK_WIDTH 6
#define BLOCK_WIDTH (1 << LOG_BLOCK_WIDTH)

static inline uint32_t intlog2(uint32_t input) {
    uint32_t output;
    frexp(input >> 1, (int *)&output);
    return output;
}

static const uint8_t bit_reverse_table[256] = {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98,
    0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 0x0C, 0x8C, 0x4C, 0xCC,
    0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2,
    0x72, 0xF2, 0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA, 0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6,
    0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE, 0x01, 0x81,
    0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1, 0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9,
    0x39, 0xB9, 0x79, 0xF9, 0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5, 0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD,
    0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD, 0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB, 0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97,
    0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF};

static inline uint32_t reverse_bits(uint32_t num, uint32_t bits) {
    uint32_t result = 0;
    unsigned int bytes = (bits + 7U) >> 3U;
    for (unsigned int i = 0; i < bytes; ++i) {
        result <<= 8U;
        result |= bit_reverse_table[num & 0xFFU];
        num >>= 8U;
    }
    result >>= (8U * bytes - bits);
    return result;
}

static inline void table_shuffle(float *real, float *imag, uint32_t log_n) {
    for (uint32_t i = 0; i < (1U << log_n); ++i) {
        uint32_t j = reverse_bits(i, log_n);
        if (j > i) {
            float tr = real[i];
            float ti = imag[i];
            real[i] = real[j];
            imag[i] = imag[j];
            real[j] = tr;
            imag[j] = ti;
        }
    }
}

static inline void cobra_shuffle(float *real, float *imag, uint32_t log_n, float *buffer_real, float *buffer_imag) {
    uint32_t num_b_bits = log_n - 2U * LOG_BLOCK_WIDTH;
    uint32_t b_size = 1U << num_b_bits;
    uint32_t shift_top = num_b_bits + LOG_BLOCK_WIDTH;

    uint32_t rev6[BLOCK_WIDTH];
    for (uint32_t i = 0; i < BLOCK_WIDTH; ++i) rev6[i] = bit_reverse_table[i] >> 2U;

    for (uint32_t b = 0; b < b_size; ++b) {
        uint32_t b_rev = reverse_bits(b, num_b_bits);
        uint32_t b_shifted = b << LOG_BLOCK_WIDTH;
        uint32_t b_rev_shifted = b_rev << LOG_BLOCK_WIDTH;

        for (uint32_t a = 0; a < BLOCK_WIDTH; ++a) {
            uint32_t a_rev = rev6[a];
            uint32_t a_shifted = a << shift_top;
            uint32_t a_rev_shifted = a_rev << LOG_BLOCK_WIDTH;
            for (uint32_t c = 0; c < BLOCK_WIDTH; ++c) {
                uint32_t idx = a_shifted | b_shifted | c;
                uint32_t buffer_idx = a_rev_shifted | c;
                buffer_real[buffer_idx] = real[idx];
                buffer_imag[buffer_idx] = imag[idx];
            }
        }

        for (uint32_t c = 0; c < BLOCK_WIDTH; ++c) {
            uint32_t c_rev = rev6[c];
            uint32_t c_rev_shifted = c_rev << shift_top;
            for (uint32_t a_rev = 0; a_rev < BLOCK_WIDTH; ++a_rev) {
                uint32_t a = rev6[a_rev];
                bool less = (a < c_rev) || (a == c_rev && b < b_rev) || (a == c_rev && b == b_rev && a_rev < c);
                if (less) {
                    uint32_t v_idx = c_rev_shifted | b_rev_shifted | a_rev;
                    uint32_t b_idx = (a_rev << LOG_BLOCK_WIDTH) | c;
                    float tr = real[v_idx];
                    float ti = imag[v_idx];
                    real[v_idx] = buffer_real[b_idx];
                    imag[v_idx] = buffer_imag[b_idx];
                    buffer_real[b_idx] = tr;
                    buffer_imag[b_idx] = ti;
                }
            }
        }

        for (uint32_t a = 0; a < BLOCK_WIDTH; ++a) {
            uint32_t a_rev = rev6[a];
            uint32_t a_shifted = a << shift_top;
            for (uint32_t c = 0; c < BLOCK_WIDTH; ++c) {
                uint32_t c_rev = rev6[c];
                bool less = (a < c_rev) || (a == c_rev && b < b_rev) || (a == c_rev && b == b_rev && a_rev < c);
                if (less) {
                    uint32_t v_idx = a_shifted | b_shifted | c;
                    uint32_t b_idx = (a_rev << LOG_BLOCK_WIDTH) | c;
                    float tr = real[v_idx];
                    float ti = imag[v_idx];
                    real[v_idx] = buffer_real[b_idx];
                    imag[v_idx] = buffer_imag[b_idx];
                    buffer_real[b_idx] = tr;
                    buffer_imag[b_idx] = ti;
                }
            }
        }
    }
}

static inline void bit_reverse_permutation(float *real, float *imag, uint32_t n, float *buffer_real, float *buffer_imag) {
    uint32_t order = intlog2(n);
    if (order <= 2U * LOG_BLOCK_WIDTH) {
        table_shuffle(real, imag, order);
    } else {
        cobra_shuffle(real, imag, order, buffer_real, buffer_imag);
    }
}

#if VECF_LEN == 16 || VECF_LEN == 8 || VECF_LEN == 4 || VECF_LEN == 2
#    define HAS_INTRA_VEC_PASS
#    if VECF_LEN == 16
#        define W16_C1 0.9238795325112867f
#        define W16_S1 0.3826834323650898f
static const VECF_INT permutations[4] __attribute__((aligned(64))) = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                                      {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15},
                                                                      {0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15},
                                                                      {0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15}};
static const VECF_INT inv_permutations[4] __attribute__((aligned(64))) = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                                          {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15},
                                                                          {0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15},
                                                                          {0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15}};
static const VECF real_twiddles[4]
    __attribute__((aligned(64))) = {{1.0f, W16_C1, (float)M_SQRT1_2, W16_S1, 0.0f, -W16_S1, -(float)M_SQRT1_2, -W16_C1, 1.0f, W16_C1, (float)M_SQRT1_2, W16_S1,
                                     0.0f, -W16_S1, -(float)M_SQRT1_2, -W16_C1},
                                    {1.0f, (float)M_SQRT1_2, 0.0f, -(float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2, 0.0f, -(float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2,
                                     0.0f, -(float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2, 0.0f, -(float)M_SQRT1_2},
                                    {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f},
                                    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}};
static const VECF imag_twiddles[4]
    __attribute__((aligned(64))) = {{0.0f, W16_S1, (float)M_SQRT1_2, W16_C1, 1.0f, W16_C1, (float)M_SQRT1_2, W16_S1, 0.0f, W16_S1, (float)M_SQRT1_2, W16_C1,
                                     1.0f, W16_C1, (float)M_SQRT1_2, W16_S1},
                                    {0.0f, (float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2, 0.0f, (float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2, 0.0f, (float)M_SQRT1_2,
                                     1.0f, (float)M_SQRT1_2, 0.0f, (float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2},
                                    {0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f},
                                    {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}};
static inline void nanofft_vec_shuffle(VECF *a, VECF *b) {
    VECF tmp = __builtin_shuffle(*a, *b, (VECF_INT){0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23});
    *b = __builtin_shuffle(*a, *b, (VECF_INT){8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31});
    *a = tmp;
}
#    elif VECF_LEN == 8
static const VECF_INT permutations[3] __attribute__((aligned(64))) = {{0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 4, 5, 2, 3, 6, 7}, {0, 2, 4, 6, 1, 3, 5, 7}};
static const VECF_INT inv_permutations[3] __attribute__((aligned(64))) = {{0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 4, 5, 2, 3, 6, 7}, {0, 4, 1, 5, 2, 6, 3, 7}};
static const VECF real_twiddles[3]
    __attribute__((aligned(64))) = {{1.0f, (float)M_SQRT1_2, 0.0f, -(float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2, 0.0f, -(float)M_SQRT1_2},
                                    {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f},
                                    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}};
static const VECF imag_twiddles[3]
    __attribute__((aligned(64))) = {{0.0f, (float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2, 0.0f, (float)M_SQRT1_2, 1.0f, (float)M_SQRT1_2},
                                    {0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f},
                                    {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}};
static inline void nanofft_vec_shuffle(VECF *a, VECF *b) {
    VECF tmp = __builtin_shuffle(*a, *b, (VECF_INT){0, 1, 2, 3, 8, 9, 10, 11});
    *b = __builtin_shuffle(*a, *b, (VECF_INT){4, 5, 6, 7, 12, 13, 14, 15});
    *a = tmp;
}
#    elif VECF_LEN == 4
static const VECF_INT permutations[2] __attribute__((aligned(64))) = {{0, 1, 2, 3}, {0, 2, 1, 3}};
static const VECF_INT inv_permutations[2] __attribute__((aligned(64))) = {{0, 1, 2, 3}, {0, 2, 1, 3}};
static const VECF real_twiddles[2] __attribute__((aligned(64))) = {{1.0f, 0.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 1.0f, 1.0f}};
static const VECF imag_twiddles[2] __attribute__((aligned(64))) = {{0.0f, 1.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 0.0f, 0.0f}};
static inline void nanofft_vec_shuffle(VECF *a, VECF *b) {
    VECF tmp = __builtin_shuffle(*a, *b, (VECF_INT){0, 1, 4, 5});
    *b = __builtin_shuffle(*a, *b, (VECF_INT){2, 3, 6, 7});
    *a = tmp;
}
#    elif VECF_LEN == 2
static const VECF_INT permutations[1] __attribute__((aligned(64))) = {{0, 1}};
static const VECF_INT inv_permutations[1] __attribute__((aligned(64))) = {{0, 1}};
static const VECF real_twiddles[1] __attribute__((aligned(64))) = {{1.0f, 1.0f}};
static const VECF imag_twiddles[1] __attribute__((aligned(64))) = {{0.0f, 0.0f}};
static inline void nanofft_vec_shuffle(VECF *a, VECF *b) {
    VECF tmp = __builtin_shuffle(*a, *b, (VECF_INT){0, 2});
    *b = __builtin_shuffle(*a, *b, (VECF_INT){1, 3});
    *a = tmp;
}
#    endif
static inline void nanofft_vec_perm(VECF *a, VECF *b, uint32_t idx) {
    *a = __builtin_shuffle(*a, permutations[idx]);
    *b = __builtin_shuffle(*b, permutations[idx]);
}
static inline void nanofft_vec_inv_perm(VECF *a, VECF *b, uint32_t idx) {
    *a = __builtin_shuffle(*a, inv_permutations[idx]);
    *b = __builtin_shuffle(*b, inv_permutations[idx]);
}
#endif

static void generate_fft_buffer(uint32_t n, float *real_buffer, float *imag_buffer) {
    uint32_t shift = 0;
    for (uint32_t step = n; step > 1U; step >>= 1U) {
        uint32_t half_step = step >> 1U;
        for (uint32_t j = 0; j < half_step; ++j) {
            float angle = -(float)j / (float)step;
            real_buffer[shift + j] = cos2pif_tls(angle);
            imag_buffer[shift + j] = -sin2pif_tls(angle);
        }
        shift += half_step;
    }
}

static void sande_tukey_in_place(float *real, float *imag, const float *twiddle_real, const float *twiddle_imag, uint32_t n) {
    uint32_t shift = 0;
    for (uint32_t step = n; step > VECF_LEN; step >>= 1U) {
        uint32_t half_step = step >> 1U;
        for (uint32_t i = 0; i < n; i += step) {
            for (uint32_t j = 0; j < half_step; j += VECF_LEN) {
                VECF r_even = LOAD_VEC(&real[i + j]);
                VECF i_even = LOAD_VEC(&imag[i + j]);
                VECF r_odd = LOAD_VEC(&real[i + j + half_step]);
                VECF i_odd = LOAD_VEC(&imag[i + j + half_step]);

                STORE_VEC(&real[i + j], r_even + r_odd);
                STORE_VEC(&imag[i + j], i_even + i_odd);

                VECF r_tmp = r_even - r_odd;
                VECF i_tmp = i_even - i_odd;
                VECF b_real = LOAD_VEC(&twiddle_real[shift + j]);
                VECF b_imag = LOAD_VEC(&twiddle_imag[shift + j]);

                STORE_VEC(&real[i + j + half_step], r_tmp * b_real - i_tmp * b_imag);
                STORE_VEC(&imag[i + j + half_step], r_tmp * b_imag + i_tmp * b_real);
            }
        }
        shift += half_step;
    }

#ifdef HAS_INTRA_VEC_PASS
    for (uint32_t i = int_max(0, (int32_t)intlog2(VECF_LEN) - (int32_t)intlog2(n)); i < (uint32_t)intlog2(VECF_LEN); ++i) {
        for (uint32_t j = 0; j < n; j += VECF_LEN * 2U) {
            VECF r_even = LOAD_VEC(&real[j]);
            VECF r_odd = LOAD_VEC(&real[j + VECF_LEN]);
            nanofft_vec_perm(&r_even, &r_odd, i);
            nanofft_vec_shuffle(&r_even, &r_odd);

            VECF i_even = LOAD_VEC(&imag[j]);
            VECF i_odd = LOAD_VEC(&imag[j + VECF_LEN]);
            nanofft_vec_perm(&i_even, &i_odd, i);
            nanofft_vec_shuffle(&i_even, &i_odd);

            VECF r_tmp = r_even - r_odd;
            VECF i_tmp = i_even - i_odd;
            r_even = r_even + r_odd;
            i_even = i_even + i_odd;
            r_odd = r_tmp * real_twiddles[i] - i_tmp * imag_twiddles[i];
            i_odd = r_tmp * imag_twiddles[i] + i_tmp * real_twiddles[i];

            nanofft_vec_shuffle(&r_even, &r_odd);
            nanofft_vec_inv_perm(&r_even, &r_odd, i);
            STORE_VEC(&real[j], r_even);
            STORE_VEC(&real[j + VECF_LEN], r_odd);

            nanofft_vec_shuffle(&i_even, &i_odd);
            nanofft_vec_inv_perm(&i_even, &i_odd, i);
            STORE_VEC(&imag[j], i_even);
            STORE_VEC(&imag[j + VECF_LEN], i_odd);
        }
    }
#else
    for (uint32_t step = uint_min(VECF_LEN, n); step > 1U; step >>= 1U) {
        uint32_t half_step = step >> 1U;
        for (uint32_t i = 0; i < n; i += step) {
            for (uint32_t j = 0; j < half_step; ++j) {
                float re = real[i + j];
                float ie = imag[i + j];
                float ro = real[i + j + half_step];
                float io = imag[i + j + half_step];
                real[i + j] = re + ro;
                imag[i + j] = ie + io;
                float rt = re - ro;
                float it = ie - io;
                real[i + j + half_step] = rt * twiddle_real[shift + j] - it * twiddle_imag[shift + j];
                imag[i + j + half_step] = rt * twiddle_imag[shift + j] + it * twiddle_real[shift + j];
            }
        }
        shift += half_step;
    }
#endif
}

static nufft1_nanofft_plan *nufft1_nanofft_make_plan(uint32_t n, float *twiddle_real, float *twiddle_imag, bool precomputed_twiddles) {
    if (!is_power_of_two(n) || !twiddle_real || !twiddle_imag) return NULL;
    nufft1_nanofft_plan *plan = (nufft1_nanofft_plan *)calloc(1, sizeof(*plan));
    if (!plan) return NULL;
    plan->n = n;
    plan->twiddle_real = twiddle_real;
    plan->twiddle_imag = twiddle_imag;
    if (!precomputed_twiddles) generate_fft_buffer(n, plan->twiddle_real, plan->twiddle_imag);
    return plan;
}

static void nufft1_nanofft_destroy_plan(nufft1_nanofft_plan *plan) {
    if (!plan) return;
    free(plan);
}

static void nufft1_nanofft_execute(const nufft1_nanofft_plan *plan, float *real, float *imag, float *cobra_real, float *cobra_imag) {
    sande_tukey_in_place(real, imag, plan->twiddle_real, plan->twiddle_imag, plan->n);
    bit_reverse_permutation(real, imag, plan->n, cobra_real, cobra_imag);
}

/* -------------------------------------------------------------------------
 * PSWF kernel, precompute, and execution
 * ------------------------------------------------------------------------- */

static const float pswf21_coeffs[PSWF21_DEGREE + 1] = {
    1.0000001350000001f,  0.38030810120000003f, -2.099454572f,        -1.930650019f,       1.5911306080000001f,
    0.41027796129999999f, 5.6755215970000004f,  -7.1979701699999996f, 2.1704596899999999f,
};

static const float pswf43_coeffs[PSWF43_DEGREE + 1] = {
    1.000000159e+00f, 3.814843807e-01f, -1.705772933e+00f, -1.574809829e+00f, 8.972535502e-01f,
    3.081707523e-01f, 3.170650492e+00f, -3.077828580e+00f, 6.059434412e-01f,
};

struct nufft1_plan {
    int Mpoints;
    int N;
    int Nout;
    int Nfft;
    int spread_stride;
    int spread_pad;
    int fft_alloc_len;
    int output_shift;
    int num_factors;
    nufft1_mode mode;
    float alpha;
    float *deconv;
    nufft1_nanofft_plan *nanofft_p;
};

struct nufft1_workspace {
    const nufft1_plan *plan;
    int capacity_mpoints;
    int capacity_factors;
    int capacity_spread_stride;
    int active_mpoints;
    int *spread_base_idx;
    float *spread_weight;
    float *shift_r;
    float *shift_i;
    float *fft_real;
    float *fft_imag;
    float *cobra_real;
    float *cobra_imag;
};

static inline int positive_mod_int(int value, int modulus) {
    int result = value % modulus;
    return result < 0 ? result + modulus : result;
}

static inline int round_up_int(int value, int multiple) { return ((value + multiple - 1) / multiple) * multiple; }

static void get_lege_roots(int n, float *roots, float *weights) {
    int half_n = (n + 1) / 2;
    for (int i = 1; i <= half_n; ++i) {
        float x = cos2pif_tls(0.5f * (float)(i - 0.25) / (float)(n + 0.5));
        float p1 = 0.0f;
        float p2 = 0.0f;
        float dp = 0.0f;
        float step;
        int iter = 0;
        do {
            p1 = 1.0f;
            p2 = 0.0f;
            for (int j = 1; j <= n; ++j) {
                float p0 = (((float)(2 * j - 1) * x * p1) - ((float)(j - 1) * p2)) / (float)j;
                p2 = p1;
                p1 = p0;
            }
            dp = (float)n * (x * p1 - p2) / (x * x - 1.0f);
            step = p1 / dp;
            x -= step;
        } while (fabsf(step) > 1.0e-7f && ++iter < 64);

        roots[i - 1] = x;
        weights[i - 1] = 2.0f / ((1.0f - x * x) * dp * dp);
        roots[n - i] = -x;
        weights[n - i] = weights[i - 1];
    }
}

static inline const float *pswf_coeffs(nufft1_mode mode, int *degree, float *c) {
    if (mode == NUFFT1_PSWF43) {
        *degree = PSWF43_DEGREE;
        *c = 15.658f;
        return pswf43_coeffs;
    }
    *degree = PSWF21_DEGREE;
    *c = (float)(M_PI * PSWF_W * 0.75 - 0.05);
    return pswf21_coeffs;
}

static inline float pswf0(float z, nufft1_mode mode) {
    int degree = 0;
    float c = 0.0f;
    const float *coeffs = pswf_coeffs(mode, &degree, &c);
    float z2 = z * z;
    float poly = coeffs[degree];
    for (int i = degree - 1; i >= 0; --i) poly = fmaf(poly, z2, coeffs[i]);
    return expf_remez_tls(-0.5f * c * z2) * poly;
}

static void pswf0_batch(const float *__restrict z_arr, float *__restrict out_arr, int n, nufft1_mode mode) {
    int degree = 0;
    float c = 0.0f;
    const float *coeffs = pswf_coeffs(mode, &degree, &c);
    int n_vec = n - (n % VECF_LEN);
    float neg_c_half = -0.5f * c;

    for (int i = 0; i < n_vec; i += VECF_LEN) {
        VECF v_z = LOAD_VEC(&z_arr[i]);
        VECF v_z2 = v_z * v_z;
        VECF v_poly = (VECF){} + coeffs[degree];
        for (int j = degree - 1; j >= 0; --j) v_poly = v_poly * v_z2 + ((VECF){} + coeffs[j]);
        VECF v_result = vecf_exp_tls(v_z2 * ((VECF){} + neg_c_half)) * v_poly;
        STORE_VEC(&out_arr[i], v_result);
    }

    for (int i = n_vec; i < n; ++i) out_arr[i] = pswf0(z_arr[i], mode);
}

nufft1_external_sizes nufft1_external_sizes_for_plan(int N, nufft1_mode mode) {
    nufft1_external_sizes sizes = {0};
    if (N <= 0 || (mode != NUFFT1_PSWF21 && mode != NUFFT1_PSWF43)) return sizes;
    int spread_stride = round_up_int(PSWF_W + VECF_LEN - 1, VECF_LEN);
    int nfft = next_fast_len(2 * N);
    sizes.twiddle_len = (size_t)nfft;
    sizes.fft_len = (size_t)nfft + (size_t)spread_stride;
    sizes.cobra_len = (size_t)BLOCK_WIDTH * (size_t)BLOCK_WIDTH;
    return sizes;
}

static int initialize_deconv(nufft1_plan *plan) {
    int p = plan->mode == NUFFT1_PSWF43 ? (4 * PSWF_W + 16) : (int)(1.5 * PSWF_W + 2);
    int num_gl_nodes = 2 * p;
    float *gl_nodes = alloc_aligned_float((size_t)num_gl_nodes);
    float *gl_weights = alloc_aligned_float((size_t)num_gl_nodes);
    float *precomp_vals = alloc_aligned_float((size_t)num_gl_nodes);
    if (!gl_nodes || !gl_weights || !precomp_vals) {
        free(precomp_vals);
        free(gl_nodes);
        free(gl_weights);
        return NUFFT1_UTIL_ERR_ALLOC;
    }

    get_lege_roots(num_gl_nodes, gl_nodes, gl_weights);
    pswf0_batch(gl_nodes, precomp_vals, p, plan->mode);
    for (int j = 0; j < p; ++j) precomp_vals[j] *= gl_weights[j];

    for (int k = 0; k < plan->Nout; ++k) {
        float xi_k = plan->alpha * (float)(k - plan->output_shift);
        float phi_hat_half = 0.0f;
        for (int j = 0; j < p; ++j) phi_hat_half += precomp_vals[j] * cos2pif_tls(xi_k * gl_nodes[j]);
        plan->deconv[k] = 1.0f / ((float)PSWF_W * phi_hat_half);
    }

    free(precomp_vals);
    free(gl_nodes);
    free(gl_weights);
    return NUFFT1_UTIL_OK;
}

static nufft1_plan *nufft1_initialize_shared_impl(int Mpoints, int N, int freq_factor, nufft1_mode mode, float *twiddle_real, float *twiddle_imag,
                                                  bool precomputed_twiddles) {
    if (Mpoints <= 0 || N <= 0 || freq_factor <= 0) return NULL;
    if (mode != NUFFT1_PSWF21 && mode != NUFFT1_PSWF43) return NULL;
    if (mode == NUFFT1_PSWF43 && (N % 4) != 0) return NULL;
    if (!twiddle_real || !twiddle_imag) return NULL;

    nufft1_plan *plan = (nufft1_plan *)calloc(1, sizeof(*plan));
    if (!plan) return NULL;

    plan->Mpoints = Mpoints;
    plan->N = N;
    plan->Nout = mode == NUFFT1_PSWF43 ? N + (N >> 1) : N;
    plan->Nfft = next_fast_len(2 * N);
    plan->output_shift = mode == NUFFT1_PSWF43 ? 3 * N / 4 : N / 2;
    plan->spread_stride = round_up_int(PSWF_W + VECF_LEN - 1, VECF_LEN);
    plan->spread_pad = plan->spread_stride;
    plan->fft_alloc_len = plan->Nfft + plan->spread_pad;
    plan->num_factors = freq_factor;
    plan->mode = mode;
    plan->alpha = (float)PSWF_W * 0.5f / (float)plan->Nfft;

    plan->deconv = alloc_aligned_float((size_t)plan->Nout);
    plan->nanofft_p = nufft1_nanofft_make_plan((uint32_t)plan->Nfft, twiddle_real, twiddle_imag, precomputed_twiddles);
    if (!plan->deconv || !plan->nanofft_p || initialize_deconv(plan) != NUFFT1_UTIL_OK) {
        nufft1_free_plan(plan);
        return NULL;
    }
    return plan;
}

nufft1_plan *nufft1_initialize_shared(int Mpoints, int N, int freq_factor, nufft1_mode mode, float *twiddle_real, float *twiddle_imag) {
    return nufft1_initialize_shared_impl(Mpoints, N, freq_factor, mode, twiddle_real, twiddle_imag, false);
}

nufft1_plan *nufft1_initialize_shared_precomputed(int Mpoints, int N, int freq_factor, nufft1_mode mode, float *twiddle_real, float *twiddle_imag) {
    return nufft1_initialize_shared_impl(Mpoints, N, freq_factor, mode, twiddle_real, twiddle_imag, true);
}

nufft1_workspace *nufft1_create_workspace(const nufft1_plan *plan, float *fft_real, float *fft_imag, float *cobra_real, float *cobra_imag) {
    if (!plan || !fft_real || !fft_imag || !cobra_real || !cobra_imag) return NULL;
    nufft1_workspace *workspace = (nufft1_workspace *)calloc(1, sizeof(*workspace));
    if (!workspace) return NULL;
    workspace->plan = plan;
    workspace->capacity_mpoints = plan->Mpoints;
    workspace->capacity_factors = plan->num_factors;
    workspace->capacity_spread_stride = plan->spread_stride;
    workspace->fft_real = fft_real;
    workspace->fft_imag = fft_imag;
    workspace->cobra_real = cobra_real;
    workspace->cobra_imag = cobra_imag;
    workspace->spread_base_idx = (int *)malloc((size_t)plan->Mpoints * (size_t)plan->num_factors * sizeof(int));
    workspace->spread_weight = alloc_aligned_float((size_t)plan->Mpoints * (size_t)plan->spread_stride * (size_t)plan->num_factors);
    workspace->shift_r = alloc_aligned_float((size_t)plan->Mpoints * (size_t)plan->num_factors);
    workspace->shift_i = alloc_aligned_float((size_t)plan->Mpoints * (size_t)plan->num_factors);
    if (!workspace->spread_base_idx || !workspace->spread_weight || !workspace->shift_r || !workspace->shift_i) {
        nufft1_free_workspace(workspace);
        return NULL;
    }
    return workspace;
}

int nufft1_workspace_set_plan(nufft1_workspace *workspace, const nufft1_plan *plan) {
    if (!workspace || !plan) return NUFFT1_UTIL_ERR_ARGUMENT;
    if (plan->Mpoints > workspace->capacity_mpoints || plan->num_factors > workspace->capacity_factors ||
        plan->spread_stride > workspace->capacity_spread_stride) {
        return NUFFT1_UTIL_ERR_ARGUMENT;
    }
    workspace->plan = plan;
    workspace->active_mpoints = 0;
    return NUFFT1_UTIL_OK;
}

void nufft1_precompute(nufft1_workspace *workspace, const double *x, int Mpoints, double df) {
    if (!workspace || !workspace->plan || !x || Mpoints <= 0 || Mpoints > workspace->plan->Mpoints) return;
    const nufft1_plan *plan = workspace->plan;
    workspace->active_mpoints = Mpoints;

    float z_buf[64] __attribute__((aligned(64)));
    float kval_buf[64] __attribute__((aligned(64)));

    for (int f_idx = 0; f_idx < plan->num_factors; ++f_idx) {
        int k_factor = f_idx + 1;
        size_t offset = (size_t)f_idx * (size_t)plan->Mpoints;
        size_t offset_spread = offset * (size_t)plan->spread_stride;

        float *current_shift_r = workspace->shift_r + offset;
        float *current_shift_i = workspace->shift_i + offset;
        int *current_spread_base_idx = workspace->spread_base_idx + offset;
        float *current_spread_weight = workspace->spread_weight + offset_spread;

        double out_shift = (double)plan->output_shift;
        double n_fft = (double)plan->Nfft;

        for (int m = 0; m < Mpoints; ++m) {
            double xm = x[m] * df * (double)k_factor;
            double phase = xm * out_shift;
            double phase_frac = phase - (double)((int)phase);
            current_shift_r[m] = cos2pif_tls((float)phase_frac);
            current_shift_i[m] = sin2pif_tls((float)phase_frac);

            double x_n = xm * n_fft;
            int m_left = (int)ceil(x_n - (double)PSWF_W * 0.5);
            for (int l = 0; l < PSWF_W; ++l) {
                int idx = m_left + l;
                double dist = (double)idx - x_n;
                z_buf[l] = (float)((2.0 * dist) / (double)PSWF_W);
            }

            pswf0_batch(z_buf, kval_buf, PSWF_W, plan->mode);

            int first = positive_mod_int(m_left, plan->Nfft);
            int base = first - (first % VECF_LEN);
            int lane_offset = first - base;
            size_t w_base = (size_t)m * (size_t)plan->spread_stride;
            current_spread_base_idx[m] = base;
            for (int l = 0; l < plan->spread_stride; ++l) current_spread_weight[w_base + (size_t)l] = 0.0f;
            for (int l = 0; l < PSWF_W; ++l) current_spread_weight[w_base + (size_t)lane_offset + (size_t)l] = kval_buf[l];
        }
    }
}

void nufft1_execute(const nufft1_workspace *workspace, const float *y_real, const float *y_imag, float *out_real, float *out_imag, int freq_factor) {
    if (!workspace || !workspace->plan || !y_real || !y_imag || !out_real || !out_imag) return;
    const nufft1_plan *plan = workspace->plan;
    int f_idx = freq_factor - 1;
    if (f_idx < 0 || f_idx >= plan->num_factors || workspace->active_mpoints <= 0) return;

    size_t offset = (size_t)f_idx * (size_t)plan->Mpoints;
    size_t offset_spread = offset * (size_t)plan->spread_stride;
    float *current_shift_r = workspace->shift_r + offset;
    float *current_shift_i = workspace->shift_i + offset;
    int *current_spread_base_idx = workspace->spread_base_idx + offset;
    float *current_spread_weight = workspace->spread_weight + offset_spread;

    memset(workspace->fft_real, 0, (size_t)plan->fft_alloc_len * sizeof(float));
    memset(workspace->fft_imag, 0, (size_t)plan->fft_alloc_len * sizeof(float));

    for (int m = 0; m < workspace->active_mpoints; ++m) {
        float yr = y_real[m];
        float yi = y_imag[m];
        float sr = current_shift_r[m];
        float si = current_shift_i[m];
        float sy_r = yr * sr - yi * si;
        float sy_i = yr * si + yi * sr;
        int base = current_spread_base_idx[m];
        size_t w_base = (size_t)m * (size_t)plan->spread_stride;

        VECF v_sy_r = (VECF){} + sy_r;
        VECF v_sy_i = (VECF){} + sy_i;
        for (int l = 0; l < plan->spread_stride; l += VECF_LEN) {
            int idx = base + l;
            VECF kval = LOAD_VEC(&current_spread_weight[w_base + (size_t)l]);
            VECF gr = LOAD_VEC(&workspace->fft_real[idx]);
            VECF gi = LOAD_VEC(&workspace->fft_imag[idx]);
            STORE_VEC(&workspace->fft_real[idx], gr + v_sy_r * kval);
            STORE_VEC(&workspace->fft_imag[idx], gi + v_sy_i * kval);
        }
    }

    for (int i = 0; i < plan->spread_pad; ++i) {
        int src = plan->Nfft + i;
        int dst = i % plan->Nfft;
        workspace->fft_real[dst] += workspace->fft_real[src];
        workspace->fft_imag[dst] += workspace->fft_imag[src];
    }

    nufft1_nanofft_execute(plan->nanofft_p, workspace->fft_real, workspace->fft_imag, workspace->cobra_real, workspace->cobra_imag);

    int split = plan->output_shift;
    int src_left = plan->Nfft - split;
    if ((split % VECF_LEN) == 0) {
        int n_vec_left = split - (split % VECF_LEN);
        int k = 0;
        for (; k < n_vec_left; k += VECF_LEN) {
            VECF scale = LOAD_VEC(&plan->deconv[k]);
            storeu_vecf(&out_real[k], LOAD_VEC(&workspace->fft_real[src_left + k]) * scale);
            storeu_vecf(&out_imag[k], LOAD_VEC(&workspace->fft_imag[src_left + k]) * scale);
        }
        for (; k < split; ++k) {
            int fft_idx = src_left + k;
            out_real[k] = workspace->fft_real[fft_idx] * plan->deconv[k];
            out_imag[k] = workspace->fft_imag[fft_idx] * plan->deconv[k];
        }

        int tail = plan->Nout - split;
        int tail_vec = tail - (tail % VECF_LEN);
        k = 0;
        for (; k < tail_vec; k += VECF_LEN) {
            VECF scale = LOAD_VEC(&plan->deconv[split + k]);
            storeu_vecf(&out_real[split + k], LOAD_VEC(&workspace->fft_real[k]) * scale);
            storeu_vecf(&out_imag[split + k], LOAD_VEC(&workspace->fft_imag[k]) * scale);
        }
        for (; k < tail; ++k) {
            out_real[split + k] = workspace->fft_real[k] * plan->deconv[split + k];
            out_imag[split + k] = workspace->fft_imag[k] * plan->deconv[split + k];
        }
    } else {
        for (int out_idx = 0; out_idx < plan->Nout; ++out_idx) {
            int k_prime = out_idx - plan->output_shift;
            int fft_idx = k_prime >= 0 ? k_prime : plan->Nfft + k_prime;
            out_real[out_idx] = workspace->fft_real[fft_idx] * plan->deconv[out_idx];
            out_imag[out_idx] = workspace->fft_imag[fft_idx] * plan->deconv[out_idx];
        }
    }
}

void nufft1_free_workspace(nufft1_workspace *workspace) {
    if (!workspace) return;
    free(workspace->spread_base_idx);
    free(workspace->spread_weight);
    free(workspace->shift_r);
    free(workspace->shift_i);
    free(workspace);
}

void nufft1_free_plan(nufft1_plan *plan) {
    if (!plan) return;
    free(plan->deconv);
    nufft1_nanofft_destroy_plan(plan->nanofft_p);
    free(plan);
}

int nufft1_plan_len(const nufft1_plan *plan) { return plan ? plan->N : 0; }

int nufft1_output_len(const nufft1_plan *plan) { return plan ? plan->Nout : 0; }

int nufft1_mpoints(const nufft1_plan *plan) { return plan ? plan->Mpoints : 0; }

/* -------------------------------------------------------------------------
 * Plan-size and twiddle-ladder helpers
 * ------------------------------------------------------------------------- */

static int bitceil_double(double x) {
    if (x <= 1.0) return 1;
    double capped = x > (double)UINT32_MAX ? (double)UINT32_MAX : x;
    return bitceil_int((unsigned int)ceil(capped));
}

double nufft1_approximate_cost(int N, int M, int block, int degree, double alpha, double beta, double gamma, int backend) {
    int block_eff = block;
    if (backend == PSWF_BACKEND_43) block_eff += block_eff >> 1;

    double N_eff = block_eff * ceil((double)N / (double)block_eff);
    double gamma_eff = gamma;
    if (degree > 0) gamma_eff *= (double)((2 * degree) + 1) / (double)((3 * degree) + 1);

    double cost = N_eff * pow((double)block, alpha);
    cost += beta * (N_eff - (double)block_eff) * (double)M / (double)block_eff;
    cost += gamma_eff * (double)block_eff;
    return cost;
}

int nufft1_optimize_plan_size(int N, int M, int degree, double alpha, double beta, double gamma, int backend) {
    const int min_block = 128;
    double start = pow((beta * (double)M / alpha), 1.0 / (alpha + 1.0));
    int block = bitceil_double(start);
    int n_cap = bitceil_double((double)N);

    if (block > n_cap) block = n_cap;
    if (block < min_block) block = min_block;

    double best = nufft1_approximate_cost(N, M, block, degree, alpha, beta, gamma, backend);
    while (block > min_block) {
        int next = block >> 1;
        if (next < min_block) next = min_block;
        double next_cost = nufft1_approximate_cost(N, M, next, degree, alpha, beta, gamma, backend);
        if (next_cost >= best) break;
        block = next;
        best = next_cost;
    }
    return block;
}

int nufft1_pswf43_plan_len_from_base(int base_len) {
    int plan_len = base_len;
    if (plan_len < 4) plan_len = 4;
    return (plan_len + 3) & ~3;
}

int nufft1_pswf43_output_len_for_plan(int plan_len) { return plan_len + (plan_len >> 1); }

int nufft1_twiddle_ladder_levels(int N, int block) {
    if (N <= 0 || block <= 0) return 1;
    size_t num_blocks = ((size_t)N + (size_t)block - 1U) / (size_t)block;
    if (num_blocks <= 1U) return 1;

    size_t max_advance = num_blocks - 1U;
    size_t stride = (size_t)MAX_TWIDDLE_REUSE;
    int levels = 1;
    while (max_advance >= stride) {
        ++levels;
        if (stride > SIZE_MAX / (size_t)MAX_TWIDDLE_REUSE) break;
        stride *= (size_t)MAX_TWIDDLE_REUSE;
    }
    return levels;
}

int nufft1_twiddle_ladder_carry_level(size_t next_block, int levels) {
    int level = 0;
    size_t stride = (size_t)MAX_TWIDDLE_REUSE;
    while (level + 1 < levels && next_block % stride == 0U) {
        ++level;
        if (stride > SIZE_MAX / (size_t)MAX_TWIDDLE_REUSE) break;
        stride *= (size_t)MAX_TWIDDLE_REUSE;
    }
    return level;
}

double nufft1_twiddle_ladder_advance(int block, int level) {
    double advance = (double)block;
    for (int i = 0; i < level; ++i) advance *= (double)MAX_TWIDDLE_REUSE;
    return advance;
}

/* -------------------------------------------------------------------------
 * Peak detection: three-point quadratic/Lagrange-form interpolation
 * ------------------------------------------------------------------------- */

typedef struct {
    float freq;
    float power;
    float cond;
} tlsf_peak;

static inline int float_is_finite(float v) {
    union {
        float f;
        uint32_t u;
    } bits = {v};
    return (bits.u & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

static inline int float_cmp(float a, float b) { return (a > b) - (a < b); }

static inline void quadratic_vertex(float x0, float y0, float x1, float y1, float x2, float y2, float *vx, float *vy) {
    float slope01 = (y1 - y0) / (x1 - x0);
    float slope12 = (y2 - y1) / (x2 - x1);
    float curvature = (slope12 - slope01) / (x2 - x0);
    if (curvature == 0.0f) {
        *vx = x1;
        *vy = y1;
        return;
    }
    float linear = slope01 - curvature * (x0 + x1);
    *vx = -linear / (2.0f * curvature);
    *vy = y0 + slope01 * (*vx - x0) + curvature * (*vx - x0) * (*vx - x1);
}

static inline float quadratic_value(float x, float x0, float y0, float x1, float y1, float x2, float y2) {
    float slope01 = (y1 - y0) / (x1 - x0);
    float slope12 = (y2 - y1) / (x2 - x1);
    float curvature = (slope12 - slope01) / (x2 - x0);
    return y0 + slope01 * (x - x0) + curvature * (x - x0) * (x - x1);
}

static void insert_peak(tlsf_peak peak, float *out_freq, float *out_power, float *out_cond, int has_cond, int *count, int max_peaks) {
    if (max_peaks <= 0) return;
    if (*count == max_peaks && float_cmp(peak.power, out_power[max_peaks - 1]) <= 0) return;

    int pos = *count;
    if (pos == max_peaks) {
        pos = max_peaks - 1;
    } else {
        ++(*count);
    }

    while (pos > 0 && float_cmp(peak.power, out_power[pos - 1]) > 0) {
        out_freq[pos] = out_freq[pos - 1];
        out_power[pos] = out_power[pos - 1];
        if (has_cond) out_cond[pos] = out_cond[pos - 1];
        --pos;
    }
    out_freq[pos] = peak.freq;
    out_power[pos] = peak.power;
    if (has_cond) out_cond[pos] = peak.cond;
}

int nufft1_get_peaks(const float *freq, const float *power, const float *cond, int n, int max_peaks, float threshold, float *out_freq, float *out_power,
                     float *out_cond, int *out_count) {
    if (!out_count || n < 0 || max_peaks < 0 || !float_is_finite(threshold)) return NUFFT1_UTIL_ERR_ARGUMENT;
    *out_count = 0;
    if (n > 0 && (!freq || !power)) return NUFFT1_UTIL_ERR_ARGUMENT;
    if (max_peaks > 0 && (!out_freq || !out_power || (cond && !out_cond))) return NUFFT1_UTIL_ERR_ARGUMENT;

    for (int i = 0; i < n; ++i) {
        if (!float_is_finite(freq[i])) return NUFFT1_UTIL_ERR_ARGUMENT;
        if (i > 0 && float_cmp(freq[i], freq[i - 1]) <= 0) return NUFFT1_UTIL_ERR_ARGUMENT;
    }
    if (n < 3 || max_peaks == 0) return NUFFT1_UTIL_OK;

    int count = 0;
    for (int idx = 1; idx < n - 1; ++idx) {
        float p0 = power[idx - 1];
        float p1 = power[idx];
        float p2 = power[idx + 1];
        if (!float_is_finite(p0) || !float_is_finite(p1) || !float_is_finite(p2)) continue;
        if (!(float_cmp(p1, p0) > 0 && float_cmp(p1, p2) > 0 && float_cmp(p1, threshold) > 0)) continue;

        tlsf_peak peak;
        quadratic_vertex(freq[idx - 1], p0, freq[idx], p1, freq[idx + 1], p2, &peak.freq, &peak.power);
        peak.cond = NAN;
        if (cond) {
            float c0 = cond[idx - 1];
            float c1 = cond[idx];
            float c2 = cond[idx + 1];
            if (float_is_finite(c0) && float_is_finite(c1) && float_is_finite(c2)) {
                peak.cond = quadratic_value(peak.freq, freq[idx - 1], c0, freq[idx], c1, freq[idx + 1], c2);
            }
        }
        insert_peak(peak, out_freq, out_power, out_cond, cond != NULL, &count, max_peaks);
    }

    *out_count = count;
    return NUFFT1_UTIL_OK;
}
