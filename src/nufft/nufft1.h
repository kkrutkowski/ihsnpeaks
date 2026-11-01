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
    /* One plan's twiddle region.  The fused-Radix-2 layout is not a suffix of
     * the next-larger plan's layout, so a plan cache must give every plan its
     * own twiddle_len-sized region (as src/metadata.h does). */
    size_t twiddle_len;
    size_t fft_len;
    /* Scratch space for the bit-reversal backend. One 64x64 tile buffer for the
     * scalar COBRA fallback, two when the build selects COBRAVO. Both real and
     * imaginary arrays must be allocated with this length. */
    size_t cobra_len;
} nufft1_external_sizes;

nufft1_external_sizes nufft1_external_sizes_for_plan(int N, nufft1_mode mode);
nufft1_plan *nufft1_initialize_shared(int Mpoints, int N, int freq_factor, nufft1_mode mode, float *twiddle_real, float *twiddle_imag);
/* twiddle_real/imag must already contain this plan's fused-Radix-2 layout.
 * They cannot be a suffix of a larger plan's twiddle buffer. */
nufft1_plan *nufft1_initialize_shared_precomputed(int Mpoints, int N, int freq_factor, nufft1_mode mode, float *twiddle_real, float *twiddle_imag);
nufft1_workspace *nufft1_create_workspace(const nufft1_plan *plan, float *fft_real, float *fft_imag, float *cobra_real, float *cobra_imag);
int nufft1_workspace_set_plan(nufft1_workspace *workspace, const nufft1_plan *plan);
void nufft1_precompute(nufft1_workspace *workspace, const double *x, int Mpoints, double df);
void nufft1_workspace_copy_precomputed(nufft1_workspace *dst, const nufft1_workspace *src);
int nufft1_workspace_get_active_mpoints(const nufft1_workspace *workspace);
void nufft1_execute(const nufft1_workspace *workspace, const float *y_real, const float *y_imag, float *out_real, float *out_imag, int freq_factor);

/* --- Deferred bit reversal -------------------------------------------------
 *
 * nufft1_execute_slots() skips the per-call bit reversal and returns the block
 * in "slot" order: the FFT positions the output window keeps, scaled by the
 * deconvolution.  Slot s holds the grid value of grid index
 * nufft1_slot_grid_index(); the ~1/3 of the FFT the window drops is never
 * stored.  Consumers then bit-reverse the *results* once per block instead of
 * the complex output per call:
 *
 *     nufft1_results_to_grid(ws, slot_r2, r2, count);             // one array
 *     nufft1_results_to_grid_pair(ws, slot_a, slot_b, a, b, n);   // two, one pass
 *     nufft1_slots_to_grid(ws, values, count);                    // in place
 *
 * The result array is zero-padded back to the FFT length, bit-reversed with the
 * normal COBRA/COBRAVO kernel and rotated into grid order, so the raw-NLL/Bayes
 * series is only ever evaluated on the requested grid points.
 *
 * The result helpers clobber the workspace FFT buffers (dead between two
 * nufft1_execute calls) and require nufft1_workspace_slots_ok(); the slot length
 * they consume is nufft1_output_len(). */
int nufft1_workspace_slots_ok(const nufft1_workspace *workspace);
void nufft1_execute_slots(const nufft1_workspace *workspace, const float *y_real, const float *y_imag, float *out_real, float *out_imag, int freq_factor);
void nufft1_results_to_grid(const nufft1_workspace *workspace, const float *slot_values, float *grid_out, int count);
void nufft1_results_to_grid_pair(const nufft1_workspace *workspace, const float *slot_a, const float *slot_b, float *grid_a, float *grid_b, int count);
void nufft1_slots_to_grid(const nufft1_workspace *workspace, float *values, int count);
void nufft1_free_workspace(nufft1_workspace *workspace);
void nufft1_free_plan(nufft1_plan *plan);

int nufft1_plan_len(const nufft1_plan *plan);
int nufft1_output_len(const nufft1_plan *plan);
int nufft1_mpoints(const nufft1_plan *plan);

double nufft1_approximate_cost(int N, int M, int block, int degree, double alpha, double beta, double gamma, int backend);
int nufft1_optimize_plan_size(int N, int M, int degree, double alpha, double beta, double gamma, int backend);
int nufft1_pswf43_plan_len_from_base(int base_len);
int nufft1_pswf43_output_len_for_plan(int plan_len);

enum { NUFFT1_UTIL_OK = 0, NUFFT1_UTIL_ERR_ARGUMENT = -1, NUFFT1_UTIL_ERR_ALLOC = -3 };

int nufft1_get_peaks(const float *freq, const float *power, const float *cond, int n, int max_peaks, float threshold, float *out_freq, float *out_power,
                     float *out_cond, int *out_count);

#ifdef __cplusplus
}
#endif

#endif /* NUFFT1_H */
