#ifndef NUFFT1_H
#define NUFFT1_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef VEC_BYTES
#    ifdef __AVX512F__
#        define VEC_BYTES 64
#    elif defined(__AVX__)
#        define VEC_BYTES 32
#    else
#        define VEC_BYTES 16
#    endif
#endif

#define NUFFT1_PSWF_WIDTH 8

typedef enum { NUFFT1_PSWF21 = 21, NUFFT1_PSWF43 = 43 } nufft1_mode;

typedef struct nufft1_plan nufft1_plan;
typedef struct nufft1_workspace nufft1_workspace;

typedef struct {
    size_t twiddle_len;
    size_t fft_len;
    size_t cobra_len;
} nufft1_external_sizes;

nufft1_external_sizes nufft1_external_sizes_for_plan(int N, nufft1_mode mode);
nufft1_plan *nufft1_initialize_shared(int Mpoints, int N, int freq_factor, nufft1_mode mode, float *twiddle_real, float *twiddle_imag);
nufft1_workspace *nufft1_create_workspace(const nufft1_plan *plan, float *fft_real, float *fft_imag, float *cobra_real, float *cobra_imag);
void nufft1_precompute(nufft1_workspace *workspace, const double *x, int Mpoints, double df);
void nufft1_execute(const nufft1_workspace *workspace, const float *y_real, const float *y_imag, float *out_real, float *out_imag, int freq_factor);
void nufft1_free_workspace(nufft1_workspace *workspace);
void nufft1_free_plan(nufft1_plan *plan);

int nufft1_plan_len(const nufft1_plan *plan);
int nufft1_output_len(const nufft1_plan *plan);
int nufft1_mpoints(const nufft1_plan *plan);

double nufft1_approximate_cost(int N, int M, int block, int degree, double alpha, double beta, double gamma, int backend);
int nufft1_optimize_plan_size(int N, int M, int degree, double alpha, double beta, double gamma, int backend);
int nufft1_pswf43_plan_len_from_base(int base_len);
int nufft1_pswf43_output_len_for_plan(int plan_len);
int nufft1_twiddle_ladder_levels(int N, int block);
int nufft1_twiddle_ladder_carry_level(size_t next_block, int levels);
double nufft1_twiddle_ladder_advance(int block, int level);

enum { NUFFT1_UTIL_OK = 0, NUFFT1_UTIL_ERR_ARGUMENT = -1, NUFFT1_UTIL_ERR_ALLOC = -3 };

int nufft1_get_peaks(const float *freq, const float *power, const float *cond, int n, int max_peaks, float threshold, float *out_freq, float *out_power,
                     float *out_cond, int *out_count);

#ifdef __cplusplus
}
#endif

#endif /* NUFFT1_H */
