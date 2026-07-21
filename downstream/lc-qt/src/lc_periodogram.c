/*
 * lc_periodogram.c — Single C translation unit bridging upstream ihsnpeaks into lc-qt.
 *
 * This is the "backend": it hosts BOTH the light-curve readout loaders (moved from
 * lc_readout.c) AND the multithreaded periodogram/spectrum generator. All upstream
 * header-only libraries with non-static definitions (fast_convert.h, qfits.h, sds.h,
 * fdist.h, convolution.h, kthread.h, ...) are compiled here exactly once, which avoids
 * duplicate-symbol link errors.
 *
 * Dispatch-ready: kept as a single bridge object with LC_API marking the intended public
 * entry points (see lc_periodogram.h); no -fvisibility=hidden is used. nufft1.c is a separate
 * object (not #included here), mirroring the ihsnpeaks release dispatch model for future
 * multi-microarchitecture builds.
 */

/* Upstream constants normally defined at the top of src/main.c. They are referenced by
 * header-only functions (process_path, print_help, print_peaks) that get compiled into this
 * translation unit, so they must be defined before the upstream headers are included. */
#ifndef DEFAULT_MEASUREMENT_SIZE
#    define DEFAULT_MEASUREMENT_SIZE 24
#endif
#ifndef IHSNPEAKS_VERSION
#    define IHSNPEAKS_VERSION "v1.1.0-preview"
#endif

#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "lc_periodogram.h"

/* kthread.h defines kt_forpool_t and pulls in params.h -> metadata.h -> scaling.h.
 * It must be included before process.h. */
#include <klib/kthread.h>
#include "process.h"

/* ===================================================================================
 * Light-curve readout (moved from lc_readout.c)
 * =================================================================================== */

#define LC_MAX_READ_BUF (64U * 1024U * 1024U) /* 64 MB hard cap */

/*
 * Allocate a minimal buffer_t suitable for read_dat / read_fits.
 * Only x, y, dy, readBuf, and len are populated; everything else is zeroed.
 */
static buffer_t *alloc_read_buffer(uint32_t maxLen, size_t readBufSize) {
    buffer_t *buf = calloc(1, sizeof(buffer_t));
    if (!buf) return NULL;

    buf->len = maxLen;
    buf->spectrum = sdsempty();
    buf->outBuf = sdsempty();

    buf->x = malloc((size_t)maxLen * sizeof(double));
    buf->y = malloc((size_t)maxLen * sizeof(float));
    buf->dy = malloc((size_t)maxLen * sizeof(float));

    if (!buf->x || !buf->y || !buf->dy) goto fail;

    if (readBufSize > 0) {
        buf->readBuf = malloc(readBufSize);
        if (!buf->readBuf) goto fail;
    }

    return buf;

fail:
    free(buf->x);
    free(buf->y);
    free(buf->dy);
    free(buf->readBuf);
    sdsfree(buf->spectrum);
    sdsfree(buf->outBuf);
    free(buf);
    return NULL;
}

static void free_read_buffer(buffer_t *buf) {
    if (!buf) return;
    free(buf->x);
    free(buf->y);
    free(buf->dy);
    free(buf->readBuf);
    sdsfree(buf->spectrum);
    sdsfree(buf->outBuf);
    free(buf);
}

/*
 * Transfer ownership of buffer arrays into lc_data_t.
 * The buffer_t shell is freed but the x/y/dy arrays live on inside lc_data_t.
 */
static int finalize(buffer_t *buf, lc_data_t *out) {
    if (buf->n == 0) {
        free_read_buffer(buf);
        return -1;
    }
    out->n = buf->n;
    out->x = buf->x;
    out->y = buf->y;
    out->dy = buf->dy;
    out->magnitude = buf->magnitude;
    out->_buf = buf; /* keep shell for detrend */

    /* Detach arrays from buf so free_read_buffer won't double-free */
    buf->x = NULL;
    buf->y = NULL;
    buf->dy = NULL;
    return 0;
}

int lc_load_dat(const char *path, lc_data_t *out) {
    memset(out, 0, sizeof(*out));

    /* Stat the file to determine allocation sizes */
    struct stat st;
    if (stat(path, &st) != 0 || st.st_size <= 0) return -1;

    size_t fileSize = (size_t)st.st_size;
    size_t readBufSize = fileSize + 1; /* +1 for null terminator */
    if (readBufSize > LC_MAX_READ_BUF) readBufSize = LC_MAX_READ_BUF;

    /*
     * Count newlines to get a tight upper-bound on measurement rows.
     * Each data row is terminated by '\n', so newline count >= row count.
     */
    uint32_t maxLen = 0;
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    {
        char chunk[8192];
        size_t rd;
        while ((rd = fread(chunk, 1, sizeof(chunk), fp)) > 0) {
            for (size_t i = 0; i < rd; i++) {
                if (chunk[i] == '\n') maxLen++;
            }
        }
    }
    fclose(fp);
    if (maxLen == 0) maxLen = 1;

    buffer_t *buf = alloc_read_buffer(maxLen, readBufSize);
    if (!buf) return -1;

    read_dat(path, buf);
    return finalize(buf, out);
}

int lc_load_fits(const char *path, lc_data_t *out) {
    memset(out, 0, sizeof(*out));

    /*
     * For FITS, peek at the table to get the row count.
     * We open the table just to read nr, then close and let read_fits
     * re-open it (keeps the wrapper simple and reuses upstream logic).
     */
    qfits_table *table = qfits_table_open(path, 1);
    if (!table) return -1;
    int nr = table->nr;
    qfits_table_close(table);

    if (nr <= 0) return -1;
    uint32_t maxLen = (uint32_t)nr;

    buffer_t *buf = alloc_read_buffer(maxLen, 0);
    if (!buf) return -1;

    read_fits(path, buf);
    return finalize(buf, out);
}

void lc_detrend(lc_data_t *data) {
    if (!data || data->n == 0 || !data->_buf) return;
    buffer_t *buf = (buffer_t *)data->_buf;

    /* Re-attach arrays so linregw_buffer can operate on them */
    buf->x = data->x;
    buf->y = data->y;
    buf->dy = data->dy;
    buf->n = data->n;

    /* Save the original time offset (linregw_buffer zeroes x[0]) */
    double t0 = data->x[0];

    /* Compute mean before detrending */
    double mean = 0.0;
    for (unsigned int i = 0; i < data->n; i++) {
        mean += (double)data->y[i];
    }
    mean /= (double)data->n;

    /* Weighted linear regression detrend (stores intercept in buf->magnitude) */
    linregw_buffer(buf);

    /* Restore original timescale */
    for (unsigned int i = 0; i < data->n; i++) {
        data->x[i] += t0;
    }

    /* Add mean value back so the plot shows physically meaningful levels */
    for (unsigned int i = 0; i < data->n; i++) {
        data->y[i] += (float)mean;
    }

    data->magnitude = buf->magnitude;

    /* Detach again */
    buf->x = NULL;
    buf->y = NULL;
    buf->dy = NULL;
}

void lc_free(lc_data_t *data) {
    if (!data) return;
    if (data->_buf) {
        free_read_buffer((buffer_t *)data->_buf);
        data->_buf = NULL;
    }
    free(data->x);
    free(data->y);
    free(data->dy);
    data->x = NULL;
    data->y = NULL;
    data->dy = NULL;
    data->n = 0;
}

/* ===================================================================================
 * Progress / cancellation (opaque handle, C11 atomics)
 * =================================================================================== */

struct lc_progress {
    _Atomic uint32_t done;
    _Atomic uint32_t total;
    _Atomic int cancel;
};

LC_API lc_progress_t *lc_progress_create(void) { return (lc_progress_t *)calloc(1, sizeof(lc_progress_t)); }

LC_API void lc_progress_destroy(lc_progress_t *p) { free(p); }

LC_API uint32_t lc_progress_done(lc_progress_t *p) { return p ? atomic_load(&p->done) : 0U; }

LC_API uint32_t lc_progress_total(lc_progress_t *p) { return p ? atomic_load(&p->total) : 0U; }

LC_API void lc_progress_request_cancel(lc_progress_t *p) {
    if (p) atomic_store(&p->cancel, 1);
}

LC_API void lc_progress_reset(lc_progress_t *p) {
    if (!p) return;
    atomic_store(&p->done, 0U);
    atomic_store(&p->total, 0U);
    atomic_store(&p->cancel, 0);
}

static void lc_progress_set_total(lc_progress_t *p, uint32_t total) {
    if (p) atomic_store(&p->total, total);
}

static void lc_progress_set_done(lc_progress_t *p, uint32_t done) {
    if (p) atomic_store(&p->done, done);
}

static void lc_progress_add_done(lc_progress_t *p, uint32_t count) {
    if (p) atomic_fetch_add(&p->done, count);
}

static int lc_progress_is_cancelled(lc_progress_t *p) { return p ? atomic_load(&p->cancel) : 0; }

/* ===================================================================================
 * IHS / AoV parallel sweeps with NuFFT-block-level progress tracking
 *
 * These replicate the upstream frequency-slice parallelism (execute_nufft_sweep_parallel,
 * execute_aov_sweep_parallel) but add atomic progress updates after each NuFFT block,
 * so the UI can display a live percentage even for the near-instantaneous IHS/AoV methods.
 * =================================================================================== */

/* Progress-aware serial IHS sweep (mirrors execute_nufft_sweep).
 * If skip_precompute is true, the caller has already precomputed (and this workspace copied it). */
static void lc_ihs_sweep_progress(buffer_t *buffer, double fmin, double fstep, uint32_t nfreq, lc_progress_t *progress, bool skip_precompute) {
    memset(buffer->power, 0, ((size_t)nfreq + 1U) * sizeof(float));
    if (!skip_precompute) {
        nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);
    }

    size_t num_blocks = ((size_t)nfreq + (size_t)buffer->activeOutputLen - 1U) / (size_t)buffer->activeOutputLen;
    for (int t = 0; t < buffer->terms; ++t) {
        double freq_factor = (double)(t + 1);
        fill_complex_input(buffer, fmin, freq_factor);
        fill_twiddle_ladder(buffer, fstep, freq_factor);
        reset_work_ladder(buffer);

        for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
            uint32_t base = (uint32_t)(block_idx * (size_t)buffer->activeOutputLen);
            uint32_t count = buffer->activeOutputLen;
            if (base + count > nfreq) count = nfreq - base;
            nufft1_execute(buffer->nufftWorkspace, buffer->workReal, buffer->workImag, buffer->blockReal, buffer->blockImag, t + 1);
            accumulate_power(buffer, base, count);
            if (block_idx + 1U < num_blocks) advance_work_ladder(buffer, block_idx + 1U);
            lc_progress_add_done(progress, 1U);
        }
    }
}

typedef struct {
    buffer_t *primary;
    parameters *params;
    double fmin, fstep;
    uint32_t nfreq;
    freq_slice_workset_t *worksets;
    lc_progress_t *progress;
} lc_ihs_slice_ctx_t;

static void lc_ihs_slice_worker(void *data, long i, int tid) {
    lc_ihs_slice_ctx_t *ctx = (lc_ihs_slice_ctx_t *)data;
    freq_slice_workset_t *work = &ctx->worksets[i];
    buffer_t *buf = ctx->params->buffers[i];
    buffer_t *primary = ctx->primary;

    work->buffer_idx = (int)i;
    if (work->slice_nfreq == 0) { work->status = 0; return; }

    double *own_x = buf->x; float *own_y = buf->y; float *own_dy = buf->dy; float *own_wy = buf->wy;
    buf->n = primary->n; buf->paddedLen = primary->paddedLen; buf->terms = primary->terms;
    buf->x = primary->x; buf->y = primary->y; buf->dy = primary->dy; buf->wy = primary->wy;

    const nufft_plan_cache_entry_t *entry = &ctx->params->nufftPlanCache[primary->activePlanIndex];
    if (nufft1_workspace_set_plan(buf->nufftWorkspace, entry->plan) != NUFFT1_UTIL_OK) {
        work->status = -1;
        buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
        return;
    }
    buf->activePlanIndex = primary->activePlanIndex;
    buf->activeGridLen = primary->activeGridLen;
    buf->activeOutputLen = primary->activeOutputLen;
    buf->activeLadderLevels = primary->activeLadderLevels;

    /* Copy precomputed spreading indices/weights from primary (hoisted nufft1_precompute) */
    nufft1_workspace_copy_precomputed(buf->nufftWorkspace, primary->nufftWorkspace);

    double slice_fmin = ctx->fmin + (double)work->freq_start * ctx->fstep;
    lc_ihs_sweep_progress(buf, slice_fmin, ctx->fstep, work->slice_nfreq, ctx->progress, true);
    work->status = 0;

    buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
}

static bool lc_execute_ihs_sweep_parallel(buffer_t *buffer, parameters *params, kt_forpool_t *pool, double fmin, double fstep, uint32_t nfreq,
                                          lc_progress_t *progress) {
    if (nfreq == 0) return true;
    if (!pool || pool->n_threads < 2 || !params->buffers || params->nbuffers < pool->n_threads) {
        /* Serial fallback: precompute internally */
        uint32_t outLen = buffer->activeOutputLen > 0 ? buffer->activeOutputLen : 1;
        uint32_t nb = (nfreq + outLen - 1U) / outLen;
        lc_progress_set_total(progress, (uint32_t)buffer->terms * nb);
        lc_ihs_sweep_progress(buffer, fmin, fstep, nfreq, progress, false);
        return true;
    }

    int nthreads = pool->n_threads;
    uint32_t nworksets = (uint32_t)nthreads;
    if (nworksets > nfreq) nworksets = nfreq;
    if (nworksets < 2U) nworksets = 2U;

    uint32_t per_thread_nfreq = (nfreq + nworksets - 1U) / nworksets;
    activate_target_nufft_plan(buffer, params, per_thread_nfreq);

    /* Hoist nufft1_precompute: compute once on primary, workers will copy (thread-safe) */
    nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);

    uint32_t outLen = buffer->activeOutputLen > 0 ? buffer->activeOutputLen : 1;

    uint32_t freqs_per_worker = (nfreq + nworksets - 2U) / (nworksets - 1U);
    freq_slice_workset_t *worksets = calloc(nworksets, sizeof(*worksets));
    if (!worksets) return false;

    /* Calculate total progress from actual worker slice sizes (accounts for overlap) */
    uint32_t total_blocks = 0;
    for (uint32_t w = 0; w < nworksets; ++w) {
        uint32_t start = w * (freqs_per_worker - 1U);
        uint32_t count = freqs_per_worker;
        if (start >= nfreq) { start = nfreq; count = 0; }
        else if (start + count > nfreq) count = nfreq - start;
        worksets[w].freq_start = start;
        worksets[w].slice_nfreq = count;
        worksets[w].status = -1;
        if (count > 0) {
            uint32_t slice_blocks = (count + outLen - 1U) / outLen;
            total_blocks += (uint32_t)buffer->terms * slice_blocks;
        }
    }
    lc_progress_set_total(progress, total_blocks);

    lc_ihs_slice_ctx_t ctx = {.primary = buffer, .params = params, .fmin = fmin, .fstep = fstep, .nfreq = nfreq, .worksets = worksets, .progress = progress};
    kt_forpool(pool, lc_ihs_slice_worker, &ctx, (long)nworksets);

    bool ok = true;
    for (uint32_t w = 0; w < nworksets; ++w)
        if (worksets[w].status != 0) ok = false;

    if (ok) {
        for (uint32_t w = 0; w < nworksets; ++w) {
            if (worksets[w].slice_nfreq == 0) continue;
            buffer_t *buf = params->buffers[w];
            uint32_t skip = (w == 0) ? 0U : 1U;
            memcpy(buffer->power + worksets[w].freq_start + skip, buf->power + skip, (size_t)(worksets[w].slice_nfreq - skip) * sizeof(float));
        }
    }
    free(worksets);
    return ok;
}

/* Progress-aware AoV trig-sums ladder (mirrors aov_compute_trig_sums_ladder). */
static void lc_aov_trig_sums_progress(buffer_t *buffer, double fmin, double fstep, size_t num_blocks, int max_factor, float epsilon, float ymean,
                                      bool yweighted, float *S, float *C, uint32_t total_count, lc_progress_t *progress) {
    uint32_t block_len = buffer->activeOutputLen;
    for (int q = 1; q <= max_factor; ++q) {
        double freq_factor = (double)q;
        aov_fill_complex_input(buffer, fmin, freq_factor, epsilon, ymean, yweighted);
        fill_twiddle_ladder(buffer, fstep, freq_factor);
        compute_work_at_block(buffer, 0);

        for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
            uint32_t base = (uint32_t)(block_idx * (size_t)block_len);
            uint32_t count = block_len;
            if (base + count > total_count) count = total_count - base;

            nufft1_execute(buffer->nufftWorkspace, buffer->workReal, buffer->workImag, buffer->blockReal, buffer->blockImag, q);
            memcpy(C + ((size_t)q * (size_t)total_count) + base, buffer->blockReal, (size_t)count * sizeof(float));
            memcpy(S + ((size_t)q * (size_t)total_count) + base, buffer->blockImag, (size_t)count * sizeof(float));

            if (block_idx + 1U < num_blocks) advance_work_ladder(buffer, block_idx + 1U);
            lc_progress_add_done(progress, 1U);
        }
    }
}

/* Progress-aware serial AoV sweep (mirrors execute_aov_sweep with store_power=true, scan_peaks=false). */
static int lc_aov_sweep_progress(buffer_t *buffer, parameters *params, double fmin, double fstep, uint32_t nfreq, const aov_reference_t *ref,
                                 lc_progress_t *progress, bool skip_precompute) {
    if (nfreq == 0) return 0;
    if (!aov_target_has_dof(buffer, params->nterms)) {
        for (uint32_t i = 0; i < nfreq; ++i) buffer->power[i] = 0.0f;
        return 0;
    }

    int degree = params->nterms;
    int max_factor = 2 * degree;
    uint32_t block_len = buffer->activeOutputLen;

    bool aov_arrays_owned = false;
    float *Sw, *Cw, *Syw, *Cyw, *power_arr;
    if (buffer->aovSw && buffer->aovArrayCap >= (size_t)nfreq) {
        Sw = buffer->aovSw; Cw = buffer->aovCw; Syw = buffer->aovSyw; Cyw = buffer->aovCyw; power_arr = buffer->aovPower;
    } else {
        aov_arrays_owned = true;
        Sw = (float *)aov_aligned_alloc((size_t)(max_factor + 1) * (size_t)nfreq, sizeof(float));
        Cw = (float *)aov_aligned_alloc((size_t)(max_factor + 1) * (size_t)nfreq, sizeof(float));
        Syw = (float *)aov_aligned_alloc((size_t)(degree + 1) * (size_t)nfreq, sizeof(float));
        Cyw = (float *)aov_aligned_alloc((size_t)(degree + 1) * (size_t)nfreq, sizeof(float));
        power_arr = (float *)aov_aligned_alloc((size_t)nfreq, sizeof(float));
        if (!Sw || !Cw || !Syw || !Cyw || !power_arr) {
            free(Sw); free(Cw); free(Syw); free(Cyw); free(power_arr);
            return -1;
        }
    }

    if (!skip_precompute && nufft1_workspace_get_active_mpoints(buffer->nufftWorkspace) != (int)buffer->n) {
        nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);
    }

    for (uint32_t i = 0; i < nfreq; ++i) {
        Sw[i] = 0.0f; Cw[i] = ref->ws; Syw[i] = 0.0f; Cyw[i] = ref->yws;
    }

    size_t num_blocks = ((size_t)nfreq + (size_t)block_len - 1U) / (size_t)block_len;
    lc_aov_trig_sums_progress(buffer, fmin, fstep, num_blocks, max_factor, params->epsilon, ref->ymean, false, Sw, Cw, nfreq, progress);
    lc_aov_trig_sums_progress(buffer, fmin, fstep, num_blocks, degree, params->epsilon, ref->ymean, true, Syw, Cyw, nfreq, progress);

    int status = 0;
    if (degree == 1) {
        aov_gls_impl(Sw, Cw, Syw, Cyw, (int)nfreq, ref->ws, ref->yws, ref->chi2_ref, power_arr);
    } else {
        status = aov_solve_periodogram_vec(Sw, Cw, Syw, Cyw, (int)nfreq, degree, ref->chi2_ref, power_arr);
    }

    if (status == 0) {
        for (uint32_t idx = 0; idx < nfreq; ++idx) buffer->power[idx] = power_arr[idx];
    }

    if (aov_arrays_owned) {
        free(Sw); free(Cw); free(Syw); free(Cyw); free(power_arr);
    }
    return status;
}

typedef struct {
    buffer_t *primary;
    parameters *params;
    double fmin, fstep;
    uint32_t nfreq;
    freq_slice_workset_t *worksets;
    lc_progress_t *progress;
    aov_reference_t ref;
} lc_aov_slice_ctx_t;

static void lc_aov_slice_worker(void *data, long i, int tid) {
    lc_aov_slice_ctx_t *ctx = (lc_aov_slice_ctx_t *)data;
    freq_slice_workset_t *work = &ctx->worksets[i];
    buffer_t *buf = ctx->params->buffers[i];
    buffer_t *primary = ctx->primary;
    parameters *params = ctx->params;

    work->buffer_idx = (int)i;
    if (work->slice_nfreq == 0) { work->status = 0; return; }

    double *own_x = buf->x; float *own_y = buf->y; float *own_dy = buf->dy; float *own_wy = buf->wy;
    buf->n = primary->n; buf->paddedLen = primary->paddedLen; buf->terms = primary->terms;
    buf->neff = primary->neff; buf->amp_neff = primary->amp_neff;
    buf->x = primary->x; buf->y = primary->y; buf->dy = primary->dy; buf->wy = primary->wy;

    const nufft_plan_cache_entry_t *entry = &params->nufftPlanCache[primary->activePlanIndex];
    if (nufft1_workspace_set_plan(buf->nufftWorkspace, entry->plan) != NUFFT1_UTIL_OK) {
        work->status = -1;
        buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
        return;
    }
    buf->activePlanIndex = primary->activePlanIndex;
    buf->activeGridLen = primary->activeGridLen;
    buf->activeOutputLen = primary->activeOutputLen;
    buf->activeLadderLevels = primary->activeLadderLevels;

    nufft1_workspace_copy_precomputed(buf->nufftWorkspace, primary->nufftWorkspace);

    double slice_fmin = ctx->fmin + (double)work->freq_start * ctx->fstep;
    int status = lc_aov_sweep_progress(buf, params, slice_fmin, ctx->fstep, work->slice_nfreq, &ctx->ref, ctx->progress, true);
    work->status = (status == 0) ? 0 : -1;

    buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
}

static bool lc_execute_aov_sweep_parallel(buffer_t *buffer, parameters *params, kt_forpool_t *pool, double fmin, double fstep, uint32_t nfreq,
                                          lc_progress_t *progress) {
    if (nfreq == 0) return true;
    if (!aov_target_has_dof(buffer, params->nterms)) return true;
    aov_reference_t ref;
    if (!aov_prepare_reference(buffer, params->epsilon, &ref)) return true;

    if (!pool || pool->n_threads < 2 || !params->buffers || params->nbuffers < pool->n_threads) {
        /* Serial fallback */
        uint32_t outLen = buffer->activeOutputLen > 0 ? buffer->activeOutputLen : 1;
        uint32_t nb = (nfreq + outLen - 1U) / outLen;
        lc_progress_set_total(progress, (uint32_t)(3 * params->nterms) * nb);
        nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);
        return lc_aov_sweep_progress(buffer, params, fmin, fstep, nfreq, &ref, progress, true) == 0;
    }

    int nthreads = pool->n_threads;
    uint32_t nworksets = (uint32_t)nthreads;
    if (nworksets > nfreq) nworksets = nfreq;
    if (nworksets < 2U) nworksets = 2U;

    uint32_t per_thread_nfreq = (nfreq + nworksets - 1U) / nworksets;
    activate_target_nufft_plan(buffer, params, per_thread_nfreq);

    /* Hoist nufft1_precompute: compute once on primary, workers will copy (thread-safe) */
    nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);

    uint32_t outLen = buffer->activeOutputLen > 0 ? buffer->activeOutputLen : 1;

    uint32_t freqs_per_worker = (nfreq + nworksets - 2U) / (nworksets - 1U);
    freq_slice_workset_t *worksets = calloc(nworksets, sizeof(*worksets));
    if (!worksets) return false;

    /* Calculate total progress from actual worker slice sizes (accounts for overlap) */
    uint32_t total_blocks = 0;
    for (uint32_t w = 0; w < nworksets; ++w) {
        uint32_t start = w * (freqs_per_worker - 1U);
        uint32_t count = freqs_per_worker;
        if (start >= nfreq) { start = nfreq; count = 0; }
        else if (start + count > nfreq) count = nfreq - start;
        worksets[w].freq_start = start;
        worksets[w].slice_nfreq = count;
        worksets[w].status = -1;
        if (count > 0) {
            uint32_t slice_blocks = (count + outLen - 1U) / outLen;
            total_blocks += (uint32_t)(3 * params->nterms) * slice_blocks;
        }
    }
    lc_progress_set_total(progress, total_blocks);

    for (uint32_t w = 0; w < nworksets; ++w) {
        if (worksets[w].slice_nfreq == 0) continue;
        buffer_t *buf = params->buffers[w];
        if (!buffer_ensure_aov_arrays(buf, (size_t)worksets[w].slice_nfreq, params->nterms)) {
            free(worksets);
            return false;
        }
    }

    lc_aov_slice_ctx_t ctx = {.primary = buffer, .params = params, .fmin = fmin, .fstep = fstep, .nfreq = nfreq,
                              .worksets = worksets, .progress = progress, .ref = ref};
    kt_forpool(pool, lc_aov_slice_worker, &ctx, (long)nworksets);

    bool ok = true;
    for (uint32_t w = 0; w < nworksets; ++w)
        if (worksets[w].status != 0) ok = false;

    if (ok) {
        for (uint32_t w = 0; w < nworksets; ++w) {
            if (worksets[w].slice_nfreq == 0) continue;
            buffer_t *buf = params->buffers[w];
            uint32_t skip = (w == 0) ? 0U : 1U;
            memcpy(buffer->power + worksets[w].freq_start + skip, buf->power + skip, (size_t)(worksets[w].slice_nfreq - skip) * sizeof(float));
        }
    }
    free(worksets);
    return ok;
}

/* ===================================================================================
 * GB/BLS direct evaluation grid (bridge-side loop for progress + cancellation)
 *
 * Reuses upstream per-frequency math (get_eval_likelihood) and per-thread scratch
 * (alloc_direct_scratch / free_direct_scratch / direct_grid_chunk_size); only the
 * thin workset orchestration is custom so we can report progress and honour cancel,
 * which upstream execute_direct_grid does not expose.
 * =================================================================================== */

typedef struct {
    buffer_t *scratch; /* per-thread scratch (length nthreads) */
    parameters *params;
    const eval_method_t *method;
    double fmin;
    double fstep;
    uint32_t nfreq;
    uint32_t chunk;
    float *out; /* NLL output, length nfreq */
    lc_progress_t *progress;
} lc_direct_ctx_t;

static void lc_direct_worker(void *data, long i, int tid) {
    lc_direct_ctx_t *ctx = (lc_direct_ctx_t *)data;
    buffer_t *s = &ctx->scratch[tid];

    uint32_t begin = (uint32_t)i * ctx->chunk;
    uint32_t end = begin + ctx->chunk;
    if (end > ctx->nfreq) end = ctx->nfreq;
    if (begin >= end) return;

    /* Cooperative cancellation: stop at chunk boundaries. */
    if (lc_progress_is_cancelled(ctx->progress)) return;

    for (uint32_t idx = begin; idx < end; ++idx) {
        double freq = ctx->fmin + (double)idx * ctx->fstep;
        float value = get_eval_likelihood(s, ctx->params, ctx->method, freq, NULL, NULL);
        ctx->out[idx] = float_is_finite_bits(value) ? value : 0.0f;
    }
    lc_progress_add_done(ctx->progress, end - begin);
}

/* Returns 0 on success, 1 if cancelled, <0 on error. */
static int lc_run_direct_grid(buffer_t *buffer, parameters *params, const eval_method_t *method, kt_forpool_t *pool, double fmin, double fstep,
                              uint32_t nfreq, float *out, lc_progress_t *progress) {
    if (nfreq == 0U) return 0;
    lc_progress_set_total(progress, nfreq);

    uint32_t chunk = direct_grid_chunk_size(buffer->n);
    uint32_t nwork = chunk > 0U ? (nfreq + chunk - 1U) / chunk : 0U;
    if (nwork == 0U) return 0;

    int nthreads = (pool && pool->n_threads > 1) ? pool->n_threads : 1;
    buffer_t *scratch = (buffer_t *)calloc((size_t)nthreads, sizeof(buffer_t));
    if (!scratch) return -1;

    int ok = 1;
    for (int t = 0; t < nthreads; ++t) {
        if (alloc_direct_scratch(&scratch[t], buffer, params, false) != 0) {
            ok = 0;
            break;
        }
    }

    if (ok) {
        lc_direct_ctx_t ctx = {.scratch = scratch,
                               .params = params,
                               .method = method,
                               .fmin = fmin,
                               .fstep = fstep,
                               .nfreq = nfreq,
                               .chunk = chunk,
                               .out = out,
                               .progress = progress};
        if (pool && nthreads > 1) {
            kt_forpool(pool, lc_direct_worker, &ctx, (long)nwork);
        } else {
            for (uint32_t w = 0; w < nwork; ++w) lc_direct_worker(&ctx, (long)w, 0);
        }
    }

    for (int t = 0; t < nthreads; ++t) {
        if (scratch[t].pidx || scratch[t].buf || scratch[t].allocated) free_direct_scratch(&scratch[t]);
    }
    free(scratch);

    if (!ok) return -1;
    if (lc_progress_is_cancelled(progress)) return 1;
    return 0;
}

/* ===================================================================================
 * Periodogram computation
 * =================================================================================== */

/* Count physical cores by reading /proc/cpuinfo unique (physical id, core id) pairs.
 * Falls back to sysconf(_SC_NPROCESSORS_ONLN) on failure. */
static int count_physical_cores(void) {
    FILE *fp = fopen("/proc/cpuinfo", "r");
    if (!fp) return (int)sysconf(_SC_NPROCESSORS_ONLN);
    char line[256];
    int phys_id = -1, core_id = -1;
    int count = 0;
    /* Simple heuristic: count unique (physical_id * 1024 + core_id) combos */
    int seen[4096];
    int nseen = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "physical id", 11) == 0) {
            phys_id = atoi(strchr(line, ':') + 1);
        } else if (strncmp(line, "core id", 7) == 0) {
            core_id = atoi(strchr(line, ':') + 1);
        } else if (line[0] == '\n' || line[0] == '\0') {
            if (phys_id >= 0 && core_id >= 0) {
                int key = phys_id * 1024 + core_id;
                int dup = 0;
                for (int j = 0; j < nseen; ++j) {
                    if (seen[j] == key) { dup = 1; break; }
                }
                if (!dup && nseen < 4096) { seen[nseen++] = key; count++; }
            }
            phys_id = -1;
            core_id = -1;
        }
    }
    /* Handle last entry if file doesn't end with blank line */
    if (phys_id >= 0 && core_id >= 0) {
        int key = phys_id * 1024 + core_id;
        int dup = 0;
        for (int j = 0; j < nseen; ++j) {
            if (seen[j] == key) { dup = 1; break; }
        }
        if (!dup && nseen < 4096) { count++; }
    }
    fclose(fp);
    return count > 0 ? count : (int)sysconf(_SC_NPROCESSORS_ONLN);
}

LC_API void lc_periodogram_result_free(lc_periodogram_result_t *r) {
    if (!r) return;
    free(r->nll);
    r->nll = NULL;
    r->nfreq = 0U;
}

LC_API int lc_compute_periodogram(const lc_data_t *data, const lc_periodogram_config_t *cfg, lc_periodogram_result_t *out, lc_progress_t *progress) {
    if (!data || data->n == 0U || !cfg || !out) return -1;
    memset(out, 0, sizeof(*out));
    lc_progress_reset(progress);

    /* --- 1. Build parameters (defaults via init_parameters, then override) --- */
    char arg0[] = "ihsnpeaks";
    char arg1[] = "lc.dat";
    char arg2[] = "24.0";
    char *fake_argv[3] = {arg0, arg1, arg2};
    parameters params = init_parameters(3, fake_argv);

    params.nterms = cfg->nterms > 0 ? cfg->nterms : 3;
    params.oversamplingFactor = cfg->oversampling > 0.0 ? (float)cfg->oversampling : 5.0f;
    params.fmin = (float)cfg->fmin; /* 0 => auto later via effective_grid_fmin */
    params.fmax = cfg->fmax > 0.0 ? (float)cfg->fmax : 24.0f;
    params.gridMode = (cfg->pswf == 21) ? NUFFT1_PSWF21 : NUFFT1_PSWF43;

    /* GB/BLS evaluation parameters: gbls[oversmoothing] and bls[oversmoothing,0.5,nbins] */
    params.gbAlpha = (float)cfg->oversmoothing;
    params.blsMinRelWidth = cfg->oversmoothing;
    params.blsMaxRelWidth = 0.5; /* hardcoded maximum for lc-qt */
    params.blsWidthCount = cfg->nbins > 0 ? cfg->nbins : 10;

    bool use_aov = false;
    bool direct_grid = false;
    switch (cfg->method) {
        case LC_SPEC_AOV:
            params.periodogramMethod = PERIODOGRAM_AOV;
            params.mode = 2;
            use_aov = true;
            break;
        case LC_SPEC_GB:
            params.periodogramMethod = PERIODOGRAM_IHS;
            params.gbEvalMode = GB_EVAL_GBLS;
            params.mode = 5;
            direct_grid = true;
            break;
        case LC_SPEC_BLS:
            params.periodogramMethod = PERIODOGRAM_IHS;
            params.gbEvalMode = GB_EVAL_BLS;
            params.mode = 5;
            direct_grid = true;
            break;
        case LC_SPEC_IHS:
        default:
            params.periodogramMethod = PERIODOGRAM_IHS;
            params.mode = 2;
            break;
    }

    /* --- 2. Sizes (replicate process_path without file scanning) --- */
    uint32_t n = data->n;
    params.maxLen = n;
    double time_span = (n > 1U) ? (data->x[n - 1] - data->x[0]) : 0.0;
    params.maxTimeSpan = time_span;
    params.maxSize = 4096U; /* readBuf unused for in-memory data; keep non-zero */
    params.maxFreqCount = estimate_frequency_count(&params, time_span);
    initialize_nufft_plan(&params); /* always: alloc_buffer requires nufftPlanCount > 0 */

    int result_code = -1;
    float *nll = NULL;
    buffer_t buf;
    memset(&buf, 0, sizeof(buf));
    int primary_allocated = 0;
    kt_forpool_t *pool = NULL;
    int workers_allocated = 0;

    /* --- 3. Primary buffer + preprocess --- */
    if (alloc_buffer(&buf, &params) != 0) goto cleanup;
    primary_allocated = 1;
    for (uint32_t i = 0; i < n; ++i) {
        buf.x[i] = data->x[i];
        buf.y[i] = data->y[i];
        buf.dy[i] = data->dy[i];
    }
    buf.n = n;
    preprocess_buffer(&buf, params.epsilon, params.mode);

    /* --- 4. Frequency grid (mirrors process_target) --- */
    const eval_method_t *eval_method = eval_method_for_params(&params);
    double fmin = effective_grid_fmin(&params, buf.x[buf.n - 1]);
    double fstep;
    if (mode_uses_direct_eval_grid(params.mode)) {
        fstep = eval_method->direct_frequency_step(buf.n, buf.x[buf.n - 1], &params);
    } else {
        fstep = 1.0 / ((double)params.nterms * (double)params.oversamplingFactor * (double)buf.x[buf.n - 1] * 0.5);
    }
    double df = 1.0 / (double)buf.x[buf.n - 1];
    uint32_t nfreq = target_frequency_count(&params, fmin, fstep);
    if (nfreq > buf.maxFreqCount) nfreq = buf.maxFreqCount;

    nll = (float *)calloc(nfreq > 0U ? (size_t)nfreq : 1U, sizeof(float));
    if (!nll) goto cleanup;

    /* --- 5. Thread pool (one shared pool per computation) --- */
    int nthreads = cfg->nthreads > 0 ? cfg->nthreads : count_physical_cores();
    if (nthreads < 1) nthreads = 1;
    if (nthreads > 1) pool = (kt_forpool_t *)kt_forpool_init(nthreads, false);

    if (direct_grid) {
        /* GB / BLS */
        int rc = lc_run_direct_grid(&buf, &params, eval_method, pool, fmin, fstep, nfreq, nll, progress);
        if (rc == 1) {
            result_code = 1; /* cancelled */
            goto cleanup;
        }
        if (rc < 0) goto cleanup;
    } else {
        /* IHS / AoV: parallel NuFFT sweep filling buf.power[], then convert to NLL */
        if (!activate_target_nufft_plan(&buf, &params, nfreq)) goto cleanup;
        params.nbuffers = nthreads;
        alloc_buffers(&params);
        workers_allocated = 1;
        for (int i = 0; i < params.nbuffers; ++i) {
            if (alloc_buffer(params.buffers[i], &params) != 0) goto cleanup;
        }

        bool ok;
        if (use_aov) {
            ok = lc_execute_aov_sweep_parallel(&buf, &params, pool, fmin, fstep, nfreq, progress);
        } else {
            ok = lc_execute_ihs_sweep_parallel(&buf, &params, pool, fmin, fstep, nfreq, progress);
        }
        if (!ok) goto cleanup;
        if (lc_progress_is_cancelled(progress)) {
            result_code = 1; /* cancelled */
            goto cleanup;
        }

        int aov_n_eff = periodogram_effective_n(&buf);
        capture_spectrum_column(&buf, &params, eval_method, use_aov, false, aov_n_eff, fmin, fstep, nfreq, params.nterms, nll, 0U);
    }

    out->nfreq = nfreq;
    out->fmin = fmin;
    out->fstep = fstep;
    out->nll = nll;
    nll = NULL; /* ownership transferred */
    result_code = 0;

cleanup:
    if (pool) kt_forpool_destroy(pool);
    if (workers_allocated) {
        for (int i = 0; i < params.nbuffers; ++i) {
            if (params.buffers[i]) {
                free_buffer(params.buffers[i]);
                free(params.buffers[i]);
            }
        }
        free(params.buffers);
        params.buffers = NULL;
    }
    if (primary_allocated) free_buffer(&buf);
    free(nll);
    free_nufft_plan_cache(&params);
    free(params.target);
    free_targets(&params.targets);
    return result_code;
}
