#ifndef CUNUFFT_H
#define CUNUFFT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { CUNUFFT_PSWF21 = 21, CUNUFFT_PSWF43 = 43 } cunufft_mode;

typedef struct cunufft_plan cunufft_plan;
typedef struct cunufft_workspace cunufft_workspace;

/* Initialize a shared GPU NUFFT plan. Allocation of twiddle buffers is handled internally. */
cunufft_plan *cunufft_initialize(int Mpoints, int N, int freq_factor, cunufft_mode mode);

/* Create a GPU NUFFT workspace with associated device memory buffers. */
cunufft_workspace *cunufft_create_workspace(const cunufft_plan *plan);

/* Perform coordinate and weight precomputations on the host and upload parameters to the GPU. */
void cunufft_precompute(cunufft_workspace *workspace, const double *x, int Mpoints, double df);

/* Execute the batched GPU NUFFT concurrently for all factors.
   Input pointers are host arrays of size Mpoints (shared across all factors);
   output pointers are host arrays of size (num_factors * Nout).
   Copies to/from device are handled internally.
   Returns the raw GPU kernel execution time in milliseconds (via CUDA events). */
float cunufft_execute_batch(const cunufft_workspace *workspace, const float *y_real, const float *y_imag, float *out_real, float *out_imag);

/* Cleanup functions. */
void cunufft_free_workspace(cunufft_workspace *workspace);
void cunufft_free_plan(cunufft_plan *plan);

#ifdef __cplusplus
}
#endif

#endif /* CUNUFFT_H */
