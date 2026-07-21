#ifndef LC_PERIODOGRAM_H
#define LC_PERIODOGRAM_H

#include <stdint.h>

#include "lc_readout.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Dispatch-ready API markers.
 *
 * LC_API marks the intended public entry points of the backend. The backend is kept as a
 * single bridge object (+ separate NuFFT object) so that a future per-microarchitecture
 * variant wrapper can rename these entry points (e.g. lc_compute_periodogram_x86_64_v3)
 * and localize the rest, mirroring the ihsnpeaks release dispatch model
 * (dispatch/build_release.py). No -fvisibility=hidden is used, so LC_API is currently a
 * documentation marker (default visibility) and all symbols link normally.
 */
#if defined(__GNUC__) && __GNUC__ >= 4
#    define LC_API __attribute__((visibility("default")))
#else
#    define LC_API
#endif

typedef enum {
    LC_SPEC_IHS = 0,
    LC_SPEC_AOV,
    LC_SPEC_GB,
    LC_SPEC_BLS
} lc_spec_method_t;

typedef struct {
    lc_spec_method_t method;
    int nterms;           /* number of harmonics (-d), default 3 */
    double oversampling;  /* oversampling factor (-o), default 5 */
    double fmin;          /* min frequency (-f); 0 = auto (2/time_span) */
    double fmax;          /* max frequency (positional), default 24 */
    int pswf;             /* NuFFT backend: 43 (pswf43) or 21 (pswf21) */
    int nthreads;         /* worker threads; 0 = auto (all online CPUs) */
    double oversmoothing; /* default 0.2 -> gbAlpha (GBLS) AND blsMinRelWidth (BLS) */
    int nbins;            /* default 10 -> blsWidthCount (BLS only) */
} lc_periodogram_config_t;

typedef struct {
    uint32_t nfreq;
    double fmin;
    double fstep;
    float *nll; /* length nfreq; release with lc_periodogram_result_free */
} lc_periodogram_result_t;

/* Opaque progress/cancellation handle (C11 atomics live inside the .c file). */
typedef struct lc_progress lc_progress_t;

LC_API lc_progress_t *lc_progress_create(void);
LC_API void lc_progress_destroy(lc_progress_t *p);
LC_API uint32_t lc_progress_done(lc_progress_t *p);
LC_API uint32_t lc_progress_total(lc_progress_t *p);
LC_API void lc_progress_request_cancel(lc_progress_t *p);
LC_API void lc_progress_reset(lc_progress_t *p);

/*
 * Compute a periodogram spectrum for `data` according to `cfg`.
 * Blocking; intended to run on a worker thread.
 * Returns 0 on success, <0 on error, 1 if cancelled via `progress`.
 * On success, fills *out (caller frees with lc_periodogram_result_free).
 */
LC_API int lc_compute_periodogram(const lc_data_t *data, const lc_periodogram_config_t *cfg, lc_periodogram_result_t *out, lc_progress_t *progress);

LC_API void lc_periodogram_result_free(lc_periodogram_result_t *r);

/*
 * Persistent computation context — caches thread pool, NuFFT plans, and worker
 * buffers between runs so that repeated periodogram computations (e.g. switching
 * between IHS/AoV/GB/BLS) only pay the sweep cost, not allocation overhead.
 */
typedef struct lc_compute_ctx lc_compute_ctx_t;

LC_API lc_compute_ctx_t *lc_compute_ctx_create(int nthreads); /* 0 = auto (physical cores) */
LC_API void lc_compute_ctx_destroy(lc_compute_ctx_t *ctx);

/* Same semantics as lc_compute_periodogram but reuses cached resources in ctx. */
LC_API int lc_compute_periodogram_ctx(lc_compute_ctx_t *ctx, const lc_data_t *data, const lc_periodogram_config_t *cfg, lc_periodogram_result_t *out,
                                      lc_progress_t *progress);

#ifdef __cplusplus
}
#endif

#endif /* LC_PERIODOGRAM_H */
