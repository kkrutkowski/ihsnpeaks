/*
 * lc_period.c — Single C translation unit bridging upstream ihsnpeaks into lc-qt.
 *
 * This is the "backend": it hosts BOTH the light-curve readout loaders (moved from
 * lc_readout.c) AND the multithreaded periodogram/spectrum generator. All upstream
 * header-only libraries with non-static definitions (fast_convert.h, qfits.h, sds.h,
 * fdist.h, convolution.h, kthread.h, ...) are compiled here exactly once, which avoids
 * duplicate-symbol link errors.
 *
 * Dispatch-ready: kept as a single bridge object with LC_API marking the intended public
 * entry points (see lc_period.h); no -fvisibility=hidden is used. nufft1.c is a separate
 * object (not #included here), mirroring the ihsnpeaks release dispatch model for future
 * multi-microarchitecture builds.
 *
 * The phased model overlay (lc_compute_phased_model) lives in lc_model.c.
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

#include "lc_period.h"

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

    /* Szego + Schur detrending (stores weighted mean in buf->magnitude) */
    detrend_buffer_szego(buf, 0.001, 3);

    /* Restore original timescale */
    for (unsigned int i = 0; i < data->n; i++) {
        data->x[i] += t0;
    }

    /* Add mean value back so the plot shows physically meaningful levels */
    for (unsigned int i = 0; i < data->n; i++) {
        data->y[i] += buf->magnitude;
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

static void lc_progress_add_total(lc_progress_t *p, uint32_t extra) {
    if (p) atomic_fetch_add(&p->total, extra);
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
static bool lc_ihs_sweep_progress(buffer_t *buffer, double fmin, double fstep, uint32_t nfreq, lc_progress_t *progress, bool skip_precompute) {
    if (!buffer_ensure_power(buffer, (size_t)nfreq + 1U)) return false;
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
    return true;
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
    work->status = lc_ihs_sweep_progress(buf, slice_fmin, ctx->fstep, work->slice_nfreq, ctx->progress, true) ? 0 : -1;

    buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
}

static bool lc_execute_ihs_sweep_parallel(buffer_t *buffer, parameters *params, kt_forpool_t *pool, double fmin, double fstep, uint32_t nfreq,
                                          lc_progress_t *progress, int max_slices) {
    if (nfreq == 0) return true;
    if (!pool || pool->n_threads < 2 || !params->buffers || params->nbuffers < pool->n_threads) {
        /* Serial fallback: precompute internally */
        uint32_t outLen = buffer->activeOutputLen > 0 ? buffer->activeOutputLen : 1;
        uint32_t nb = (nfreq + outLen - 1U) / outLen;
        lc_progress_set_total(progress, (uint32_t)buffer->terms * nb);
        return lc_ihs_sweep_progress(buffer, fmin, fstep, nfreq, progress, false);
    }

    int nthreads = pool->n_threads;
    uint32_t nworksets = (uint32_t)nthreads;
    if (max_slices > 0 && nworksets > (uint32_t)max_slices) nworksets = (uint32_t)max_slices;
    if (nworksets > nfreq) nworksets = nfreq;
    if (nworksets < 2U) nworksets = 2U;

    uint32_t per_thread_nfreq = (nfreq + nworksets - 1U) / nworksets;
    (void)per_thread_nfreq;

    /* Hoist nufft1_precompute: compute once on primary, workers will copy (thread-safe).
     * NOTE: keep the primary's full plan (activated with nfreq before this sweep) instead
     * of re-activating a per-thread plan: the smaller plan/ladder changes the numerics and
     * shifts the grid vs serial. Workers copy this plan, so parallel == serial bit-for-bit. */
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


/* Progress-aware AoV trig-sums ladder (mirrors aov_compute_trig_sums_ladder): computes
 * blocks [block_begin, block_end) with the hierarchical work state at block_begin and
 * sequential advance between blocks. base_offset is the absolute freq index of the
 * start of the output arrays (global base of the current sweep block). */
static void lc_aov_trig_sums_progress(buffer_t *buffer, double fmin, double fstep, size_t block_begin, size_t block_end, int max_factor, float epsilon,
                                      float ymean, bool yweighted, float *S, float *C, uint32_t total_count, uint32_t base_offset,
                                      lc_progress_t *progress) {
    uint32_t block_len = buffer->activeOutputLen;
    for (int q = 1; q <= max_factor; ++q) {
        double freq_factor = (double)q;
        aov_fill_complex_input(buffer, fmin, freq_factor, epsilon, ymean, yweighted);
        fill_twiddle_ladder(buffer, fstep, freq_factor);
        compute_work_at_block(buffer, block_begin);

        for (size_t block_idx = block_begin; block_idx < block_end; ++block_idx) {
            uint32_t base = (uint32_t)(block_idx * (size_t)block_len);
            uint32_t count = block_len;
            if (base + count > base_offset + total_count) count = base_offset + total_count - base;
            uint32_t rel_base = base - base_offset;

            nufft1_execute(buffer->nufftWorkspace, buffer->workReal, buffer->workImag, buffer->blockReal, buffer->blockImag, q);
            memcpy(C + ((size_t)q * (size_t)total_count) + rel_base, buffer->blockReal, (size_t)count * sizeof(float));
            memcpy(S + ((size_t)q * (size_t)total_count) + rel_base, buffer->blockImag, (size_t)count * sizeof(float));

            if (block_idx + 1U < block_end) advance_work_ladder(buffer, block_idx + 1U);
            lc_progress_add_done(progress, 1U);
        }
    }
}

/* Progress-aware serial AoV sweep (mirrors execute_aov_sweep with store_power=true, scan_peaks=false). */
static int lc_aov_sweep_progress(buffer_t *buffer, parameters *params, double fmin, double fstep, uint32_t nfreq, const aov_reference_t *ref,
                                 lc_progress_t *progress, bool skip_precompute) {
    if (nfreq == 0) return 0;
    if (!buffer_ensure_power(buffer, (size_t)nfreq + 1U)) return -1;
    if (!aov_target_has_dof(buffer, params->nterms)) {
        for (uint32_t i = 0; i < nfreq; ++i) buffer->power[i] = 0.0f;
        return 0;
    }

    int degree = params->nterms;
    int max_factor = 2 * degree;
    uint32_t block_len = buffer->activeOutputLen;
    if (block_len == 0U) return -1;

    // Block-sized arrays (stride = count <= block_len): each outputLen-sized block is
    // filled with trig sums, solved, and stored, then its arrays are reused by the
    // next block, so memory scales with one FFT-grid block instead of the full grid.
    if (!buffer_ensure_aov_arrays(buffer, (size_t)block_len, params->nterms)) return -1;
    float *Sw = buffer->aovSw, *Cw = buffer->aovCw, *Syw = buffer->aovSyw, *Cyw = buffer->aovCyw, *power_arr = buffer->aovPower;

    if (!skip_precompute && nufft1_workspace_get_active_mpoints(buffer->nufftWorkspace) != (int)buffer->n) {
        nufft1_precompute(buffer->nufftWorkspace, buffer->x, (int)buffer->n, fstep);
    }

    size_t num_blocks = ((size_t)nfreq + (size_t)block_len - 1U) / (size_t)block_len;
    for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
        uint32_t base = (uint32_t)(block_idx * (size_t)block_len);
        uint32_t count = block_len;
        if (base + count > nfreq) count = nfreq - base;

        for (uint32_t i = 0; i < count; ++i) {
            Sw[i] = 0.0f; Cw[i] = ref->ws; Syw[i] = 0.0f; Cyw[i] = ref->yws;
        }

        lc_aov_trig_sums_progress(buffer, fmin, fstep, block_idx, block_idx + 1U, max_factor, params->epsilon, ref->ymean, false, Sw, Cw, count, base,
                                  progress);
        lc_aov_trig_sums_progress(buffer, fmin, fstep, block_idx, block_idx + 1U, degree, params->epsilon, ref->ymean, true, Syw, Cyw, count, base,
                                  progress);

        float *cond_arr = buffer->aovCondition;
        if (degree == 1) {
            aov_gls_impl(Sw, Cw, Syw, Cyw, (int)count, ref->ws, ref->yws, ref->chi2_ref, power_arr, cond_arr);
        } else {
            int status = aov_solve_periodogram_vec(Sw, Cw, Syw, Cyw, (int)count, degree, ref->chi2_ref, power_arr, cond_arr);
            if (status != 0) return status;
        }

        if (params->statistic == STATISTIC_BAYES) {
            int n_eff = periodogram_effective_n(buffer);
            nll_convert_spectrum_bayes_batch(power_arr, cond_arr, power_arr, count, degree, n_eff);
        }

        for (uint32_t i = 0; i < count; ++i) buffer->power[base + i] = power_arr[i];
    }
    return 0;
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
    work->power_direct = false;
    if (work->slice_nfreq == 0) { work->status = 0; return; }

    double *own_x = buf->x; float *own_y = buf->y; float *own_dy = buf->dy; float *own_wy = buf->wy;
    float *own_power = buf->power;
    size_t own_powerCap = buf->powerCap;
    buf->n = primary->n; buf->paddedLen = primary->paddedLen; buf->terms = primary->terms;
    buf->neff = primary->neff; buf->amp_neff = primary->amp_neff;
    buf->x = primary->x; buf->y = primary->y; buf->dy = primary->dy; buf->wy = primary->wy;

    const nufft_plan_cache_entry_t *entry = &params->nufftPlanCache[primary->activePlanIndex];
    if (nufft1_workspace_set_plan(buf->nufftWorkspace, entry->plan) != NUFFT1_UTIL_OK) {
        work->status = -1;
        buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
        buf->power = own_power; buf->powerCap = own_powerCap;
        return;
    }
    buf->activePlanIndex = primary->activePlanIndex;
    buf->activeGridLen = primary->activeGridLen;
    buf->activeOutputLen = primary->activeOutputLen;
    buf->activeLadderLevels = primary->activeLadderLevels;

    nufft1_workspace_copy_precomputed(buf->nufftWorkspace, primary->nufftWorkspace);

    // Aggressive path: workers write their power slice directly into the primary's
    // full power grid at global offsets (disjoint ranges, 1-freq overlap at boundaries
    // gets the same value rewritten by both neighbors). Borrow at freq_start (NOT
    // freq_start+skip): lc_aov_sweep_progress writes all local indices [0, slice_nfreq)
    // into buf->power, so local i must map to global freq_start+i.
    bool power_borrowed = false;
    if (primary->power) {
        size_t borrowed = (size_t)work->freq_start;
        if (primary->powerCap >= borrowed + (size_t)work->slice_nfreq + 1U) {
            buf->power = primary->power + borrowed;
            buf->powerCap = primary->powerCap - borrowed;
            power_borrowed = true;
        }
    }
    work->power_direct = power_borrowed;

    double slice_fmin = ctx->fmin + (double)work->freq_start * ctx->fstep;
    int status = lc_aov_sweep_progress(buf, params, slice_fmin, ctx->fstep, work->slice_nfreq, &ctx->ref, ctx->progress, true);
    work->status = (status == 0) ? 0 : -1;

    if (power_borrowed) { buf->power = own_power; buf->powerCap = own_powerCap; }
    buf->x = own_x; buf->y = own_y; buf->dy = own_dy; buf->wy = own_wy;
}

static bool lc_execute_aov_sweep_parallel(buffer_t *buffer, parameters *params, kt_forpool_t *pool, double fmin, double fstep, uint32_t nfreq,
                                          lc_progress_t *progress, int max_slices) {
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
    if (max_slices > 0 && nworksets > (uint32_t)max_slices) nworksets = (uint32_t)max_slices;
    if (nworksets > nfreq) nworksets = nfreq;
    if (nworksets < 2U) nworksets = 2U;

    uint32_t per_thread_nfreq = (nfreq + nworksets - 1U) / nworksets;
    (void)per_thread_nfreq;

    /* Hoist nufft1_precompute: compute once on primary, workers will copy (thread-safe).
     * NOTE: keep the primary's full plan (activated with nfreq before this sweep) instead
     * of re-activating a per-thread plan: the smaller plan/ladder changes the numerics and
     * shifts the grid vs serial. Workers copy this plan, so parallel == serial bit-for-bit. */
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

    // AoV arrays are ensured block-sized inside lc_aov_sweep_progress (per-worker,
    // on demand); allocation failure surfaces via work->status.

    lc_aov_slice_ctx_t ctx = {.primary = buffer, .params = params, .fmin = fmin, .fstep = fstep, .nfreq = nfreq,
                              .worksets = worksets, .progress = progress, .ref = ref};
    kt_forpool(pool, lc_aov_slice_worker, &ctx, (long)nworksets);

    bool ok = true;
    for (uint32_t w = 0; w < nworksets; ++w)
        if (worksets[w].status != 0) ok = false;

    if (ok) {
        // Workers that wrote directly into primary->power are already in place.
        for (uint32_t w = 0; w < nworksets; ++w) {
            if (worksets[w].slice_nfreq == 0 || worksets[w].power_direct) continue;
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
 * =================================================================================== */

typedef struct {
    buffer_t *scratch;
    parameters *params;
    const eval_method_t *method;
    double fmin;
    double fstep;
    uint32_t nfreq;
    uint32_t chunk;
    float *out;
    lc_progress_t *progress;
} lc_direct_ctx_t;

static void lc_direct_worker(void *data, long i, int tid) {
    lc_direct_ctx_t *ctx = (lc_direct_ctx_t *)data;
    buffer_t *s = &ctx->scratch[tid];

    uint32_t begin = (uint32_t)i * ctx->chunk;
    uint32_t end = begin + ctx->chunk;
    if (end > ctx->nfreq) end = ctx->nfreq;
    if (begin >= end) return;

    if (lc_progress_is_cancelled(ctx->progress)) return;

    for (uint32_t idx = begin; idx < end; ++idx) {
        double freq = ctx->fmin + (double)idx * ctx->fstep;
        float value = get_eval_likelihood(s, ctx->params, ctx->method, freq, NULL, NULL);
        ctx->out[idx] = float_is_finite_bits(value) ? value : 0.0f;
    }
    lc_progress_add_done(ctx->progress, end - begin);
}

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
        lc_direct_ctx_t ctx = {.scratch = scratch, .params = params, .method = method, .fmin = fmin, .fstep = fstep,
                               .nfreq = nfreq, .chunk = chunk, .out = out, .progress = progress};
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
 * Parallel NLL conversion (power[] -> NLL)
 * =================================================================================== */

typedef struct {
    const float *power;
    float *nll;
    uint32_t count;
    int nterms;
    int aov_n_eff;
    bool use_aov;
    bool is_bayes;
    lc_progress_t *progress;
} lc_nll_workset_t;

typedef struct {
    lc_nll_workset_t *worksets;
} lc_nll_dispatch_t;

static void lc_nll_worker(void *data, long i, int tid) {
    (void)tid;
    lc_nll_dispatch_t *disp = (lc_nll_dispatch_t *)data;
    lc_nll_workset_t *ws = &disp->worksets[i];
    if (ws->use_aov) {
        if (ws->is_bayes) {
            for (uint32_t j = 0; j < ws->count; ++j) {
                float val = ws->power[j];
                ws->nll[j] = (float_is_finite_bits(val) && val > 0.0f) ? val : 0.0f;
            }
        } else {
            nll_convert_spectrum_batch(ws->power, ws->nll, ws->count, ws->nterms, ws->aov_n_eff);
            for (uint32_t j = 0; j < ws->count; ++j) {
                if (!float_is_finite_bits(ws->nll[j]) || ws->nll[j] < 0.0f) ws->nll[j] = 0.0f;
            }
        }
    } else {
        for (uint32_t j = 0; j < ws->count; ++j) {
            float magnitude = (float)correct_ihs_res((double)ws->power[j], ws->nterms);
            ws->nll[j] = (float_is_finite_bits(magnitude) && magnitude > 0.0f) ? magnitude : 0.0f;
        }
    }
    lc_progress_add_done(ws->progress, 1U);
}

static uint32_t lc_parallel_nll_convert(const float *power, float *nll, uint32_t nfreq, int nterms, int aov_n_eff, bool use_aov, bool is_bayes,
                                        kt_forpool_t *pool, lc_progress_t *progress) {
    if (nfreq == 0) return 0;
    int nthreads = (pool && pool->n_threads > 1) ? pool->n_threads : 1;
    uint32_t nwork = (uint32_t)nthreads;
    if (nwork > nfreq) nwork = nfreq;

    uint32_t chunk = (nfreq + nwork - 1U) / nwork;

    lc_nll_workset_t *worksets = (lc_nll_workset_t *)calloc(nwork, sizeof(lc_nll_workset_t));
    if (!worksets) {
        if (use_aov) {
            if (is_bayes) {
                for (uint32_t j = 0; j < nfreq; ++j) {
                    float val = power[j];
                    nll[j] = (float_is_finite_bits(val) && val > 0.0f) ? val : 0.0f;
                }
            } else {
                nll_convert_spectrum_batch(power, nll, nfreq, nterms, aov_n_eff);
                for (uint32_t j = 0; j < nfreq; ++j) {
                    if (!float_is_finite_bits(nll[j]) || nll[j] < 0.0f) nll[j] = 0.0f;
                }
            }
        } else {
            for (uint32_t j = 0; j < nfreq; ++j) {
                float magnitude = (float)correct_ihs_res((double)power[j], nterms);
                nll[j] = (float_is_finite_bits(magnitude) && magnitude > 0.0f) ? magnitude : 0.0f;
            }
        }
        return 0;
    }

    for (uint32_t w = 0; w < nwork; ++w) {
        uint32_t start = w * chunk;
        uint32_t count = (start + chunk <= nfreq) ? chunk : (nfreq - start);
        worksets[w].power = power + start;
        worksets[w].nll = nll + start;
        worksets[w].count = count;
        worksets[w].nterms = nterms;
        worksets[w].aov_n_eff = aov_n_eff;
        worksets[w].use_aov = use_aov;
        worksets[w].is_bayes = is_bayes;
        worksets[w].progress = progress;
    }

    lc_nll_dispatch_t disp = {.worksets = worksets};
    if (pool && pool->n_threads > 1) {
        kt_forpool(pool, lc_nll_worker, &disp, (long)nwork);
    } else {
        for (uint32_t w = 0; w < nwork; ++w) lc_nll_worker(&disp, (long)w, 0);
    }
    free(worksets);
    return nwork;
}

/* ===================================================================================
 * Peak detection (quadratic Lagrange, mode-0 style)
 * =================================================================================== */

static bool lc_quadratic_peak_position(double freq, double fstep, float left, float center, float right, double oversampling, double *peak_freq,
                                       float *peak_nll) {
    if (!float_is_finite_bits(left) || !float_is_finite_bits(center) || !float_is_finite_bits(right)) return false;
    *peak_freq = freq;
    *peak_nll = center;
    if (oversampling < 3.0) return true;

    double x0 = freq - fstep;
    double x1 = freq;
    double x2 = freq + fstep;
    double slope01 = ((double)center - (double)left) / (x1 - x0);
    double slope12 = ((double)right - (double)center) / (x2 - x1);
    double curvature = (slope12 - slope01) / (x2 - x0);
    if (curvature == 0.0) return true;

    double linear = slope01 - curvature * (x0 + x1);
    double vx = -linear / (2.0 * curvature);
    double vy = (double)left + slope01 * (vx - x0) + curvature * (vx - x0) * (vx - x1);
    float vyf = (float)vy;
    if (double_is_finite_bits(vx) && float_is_finite_bits(vyf) && vx >= x0 && vx <= x2) {
        *peak_freq = vx;
        *peak_nll = vyf;
    }
    return true;
}

static void lc_insert_peak(lc_peak_t *peaks, int *npeaks, int max_peaks, double freq, float nll) {
    if (*npeaks == max_peaks && nll <= peaks[max_peaks - 1].nll) return;

    int pos = *npeaks;
    if (pos == max_peaks) {
        pos = max_peaks - 1;
    } else {
        ++(*npeaks);
    }
    while (pos > 0 && peaks[pos - 1].nll < nll) {
        peaks[pos] = peaks[pos - 1];
        --pos;
    }
    peaks[pos].freq = freq;
    peaks[pos].nll = nll;
}

static void lc_detect_peaks(const float *nll, uint32_t nfreq, double fmin, double fstep, double oversampling, double threshold_cli,
                            lc_peak_t *peaks, int *npeaks) {
    *npeaks = 0;
    if (nfreq < 3) return;

    double threshold = (threshold_cli > 0.0 ? threshold_cli : 8.0) * M_LN10;

    for (uint32_t i = 1; i + 1 < nfreq; ++i) {
        float left = nll[i - 1], center = nll[i], right = nll[i + 1];
        if (center > left && center > right && (double)center > threshold) {
            double freq = fmin + (double)i * fstep;
            double peak_freq = freq;
            float peak_nll = center;
            lc_quadratic_peak_position(freq, fstep, left, center, right, oversampling, &peak_freq, &peak_nll);
            lc_insert_peak(peaks, npeaks, LC_MAX_PEAKS, peak_freq, peak_nll);
        }
    }
}

/* ===================================================================================
 * Periodogram computation
 * =================================================================================== */

static int count_physical_cores(void) {
    FILE *fp = fopen("/proc/cpuinfo", "r");
    if (!fp) return (int)sysconf(_SC_NPROCESSORS_ONLN);
    char line[256];
    int phys_id = -1, core_id = -1;
    int count = 0;
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

/* ===================================================================================
 * Persistent computation context
 * =================================================================================== */

struct lc_compute_ctx {
    int nthreads;
    int physical_cores;
    kt_forpool_t *pool;
    parameters params;
    buffer_t primary;
    int primary_allocated;
    int workers_allocated;
    uint32_t cached_maxLen;
    uint32_t cached_maxFreqCount;
    int cached_gridMode;
    int cached_nterms;
    int cached_periodogramMethod;
    float cached_fmax;
    float cached_fmin;
    float cached_oversampling;
    double cached_oversmoothing;
    int cached_nbins;
    int cached_method;
    int cached_statistic;
};

LC_API lc_compute_ctx_t *lc_compute_ctx_create(int nthreads) {
    lc_compute_ctx_t *ctx = (lc_compute_ctx_t *)calloc(1, sizeof(lc_compute_ctx_t));
    if (!ctx) return NULL;
    ctx->physical_cores = count_physical_cores();
    if (ctx->physical_cores < 1) ctx->physical_cores = 1;
    ctx->nthreads = nthreads > 0 ? nthreads : (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (ctx->nthreads < 1) ctx->nthreads = 1;
    if (ctx->nthreads > 1) ctx->pool = (kt_forpool_t *)kt_forpool_init(ctx->nthreads, false);
    return ctx;
}

LC_API void lc_compute_ctx_destroy(lc_compute_ctx_t *ctx) {
    if (!ctx) return;
    if (ctx->pool) kt_forpool_destroy(ctx->pool);
    if (ctx->workers_allocated) {
        for (int i = 0; i < ctx->params.nbuffers; ++i) {
            if (ctx->params.buffers[i]) {
                free_buffer(ctx->params.buffers[i]);
                free(ctx->params.buffers[i]);
            }
        }
        free(ctx->params.buffers);
        ctx->params.buffers = NULL;
    }
    if (ctx->primary_allocated) free_buffer(&ctx->primary);
    free_nufft_plan_cache(&ctx->params);
    free(ctx->params.target);
    free_targets(&ctx->params.targets);
    free(ctx);
}

static int ctx_ensure_resources(lc_compute_ctx_t *ctx, uint32_t n, double time_span, const lc_periodogram_config_t *cfg) {
    int gridMode = (cfg->pswf == 21) ? NUFFT1_PSWF21 : NUFFT1_PSWF43;
    int nterms = cfg->nterms > 0 ? cfg->nterms : 3;
    int periMethod = (cfg->method == LC_SPEC_AOV) ? PERIODOGRAM_AOV : PERIODOGRAM_IHS;
    float fmax = cfg->fmax > 0.0 ? (float)cfg->fmax : 24.0f;
    float fmin = (float)cfg->fmin;
    float oversampling = cfg->oversampling > 0.0 ? (float)cfg->oversampling : 5.0f;
    double oversmoothing = cfg->oversmoothing;
    int nbins = cfg->nbins > 0 ? cfg->nbins : 10;
    int statistic = (cfg->method == LC_SPEC_AOV && cfg->statistic == LC_STAT_RAW) ? STATISTIC_RAW : ((cfg->method == LC_SPEC_AOV) ? STATISTIC_BAYES : STATISTIC_RAW);

    bool need_rebuild = (ctx->cached_maxLen != n || ctx->cached_gridMode != gridMode || ctx->cached_nterms != nterms ||
                         ctx->cached_periodogramMethod != periMethod || ctx->cached_fmax != fmax ||
                         ctx->cached_fmin != fmin || ctx->cached_oversampling != oversampling ||
                         ctx->cached_oversmoothing != oversmoothing || ctx->cached_nbins != nbins ||
                         ctx->cached_method != (int)cfg->method || ctx->cached_statistic != statistic || !ctx->primary_allocated);

    if (!need_rebuild) return 0;

    if (ctx->workers_allocated) {
        for (int i = 0; i < ctx->params.nbuffers; ++i) {
            if (ctx->params.buffers[i]) {
                free_buffer(ctx->params.buffers[i]);
                free(ctx->params.buffers[i]);
            }
        }
        free(ctx->params.buffers);
        ctx->params.buffers = NULL;
        ctx->workers_allocated = 0;
    }
    if (ctx->primary_allocated) {
        free_buffer(&ctx->primary);
        ctx->primary_allocated = 0;
    }
    free_nufft_plan_cache(&ctx->params);
    free(ctx->params.target);
    ctx->params.target = NULL;
    free_targets(&ctx->params.targets);

    char arg0[] = "ihsnpeaks";
    char arg1[] = "lc.dat";
    char arg2[] = "24.0";
    char *fake_argv[3] = {arg0, arg1, arg2};
    ctx->params = init_parameters(3, fake_argv);

    ctx->params.nterms = nterms;
    ctx->params.oversamplingFactor = cfg->oversampling > 0.0 ? (float)cfg->oversampling : 5.0f;
    ctx->params.fmin = (float)cfg->fmin;
    ctx->params.fmax = cfg->fmax > 0.0 ? (float)cfg->fmax : 24.0f;
    ctx->params.gridMode = gridMode;
    ctx->params.gbAlpha = (float)cfg->oversmoothing;
    ctx->params.blsMinRelWidth = cfg->oversmoothing;
    ctx->params.blsMaxRelWidth = 0.5;
    ctx->params.blsWidthCount = cfg->nbins > 0 ? cfg->nbins : 10;
    ctx->params.periodogramMethod = periMethod;
    ctx->params.statistic = (statistic_type_t)statistic;
    ctx->params.mode = (cfg->method == LC_SPEC_GB || cfg->method == LC_SPEC_BLS) ? 5 : 0;
    if (cfg->method == LC_SPEC_GB) ctx->params.gbEvalMode = GB_EVAL_GBLS;
    if (cfg->method == LC_SPEC_BLS) ctx->params.gbEvalMode = GB_EVAL_BLS;

    ctx->params.maxLen = n;
    ctx->params.maxTimeSpan = time_span;
    ctx->params.maxSize = 4096U;
    ctx->params.maxFreqCount = estimate_frequency_count(&ctx->params, time_span);
    initialize_nufft_plan(&ctx->params);

    memset(&ctx->primary, 0, sizeof(ctx->primary));
    ctx->params.spectrum = 1;  // primary holds the full power grid (NLL conversion reads it)
    if (alloc_buffer(&ctx->primary, &ctx->params) != 0) return -1;
    ctx->primary_allocated = 1;

    int nthreads = ctx->nthreads;
    ctx->params.nbuffers = nthreads;
    alloc_buffers(&ctx->params);
    ctx->workers_allocated = 1;
    ctx->params.spectrum = 0;  // workers only need FFT-grid baseline; per-slice power grown lazily
    for (int i = 0; i < ctx->params.nbuffers; ++i) {
        if (alloc_buffer(ctx->params.buffers[i], &ctx->params) != 0) return -1;
    }
    ctx->params.spectrum = 1;  // restore: sweeps and NLL conversion run in spectrum mode

    ctx->cached_maxLen = n;
    ctx->cached_maxFreqCount = ctx->params.maxFreqCount;
    ctx->cached_gridMode = gridMode;
    ctx->cached_nterms = nterms;
    ctx->cached_periodogramMethod = periMethod;
    ctx->cached_fmax = fmax;
    ctx->cached_fmin = fmin;
    ctx->cached_oversampling = oversampling;
    ctx->cached_oversmoothing = oversmoothing;
    ctx->cached_nbins = nbins;
    ctx->cached_method = (int)cfg->method;
    ctx->cached_statistic = statistic;
    return 0;
}

LC_API int lc_compute_periodogram_ctx(lc_compute_ctx_t *ctx, const lc_data_t *data, const lc_periodogram_config_t *cfg, lc_periodogram_result_t *out,
                                      lc_progress_t *progress) {
    if (!ctx || !data || data->n == 0U || !cfg || !out) return -1;
    memset(out, 0, sizeof(*out));
    lc_progress_reset(progress);

    uint32_t n = data->n;
    double time_span = (n > 1U) ? (data->x[n - 1] - data->x[0]) : 0.0;

    if (ctx_ensure_resources(ctx, n, time_span, cfg) != 0) return -1;

    parameters *params = &ctx->params;
    buffer_t *buf = &ctx->primary;

    params->oversamplingFactor = cfg->oversampling > 0.0 ? (float)cfg->oversampling : 5.0f;
    params->fmin = (float)cfg->fmin;
    params->fmax = cfg->fmax > 0.0 ? (float)cfg->fmax : 24.0f;
    params->gbAlpha = (float)cfg->oversmoothing;
    params->blsMinRelWidth = cfg->oversmoothing;
    params->blsWidthCount = cfg->nbins > 0 ? cfg->nbins : 10;
    params->statistic = (cfg->method == LC_SPEC_AOV && cfg->statistic == LC_STAT_RAW) ? STATISTIC_RAW : ((cfg->method == LC_SPEC_AOV) ? STATISTIC_BAYES : STATISTIC_RAW);
    params->mode = (cfg->method == LC_SPEC_GB || cfg->method == LC_SPEC_BLS) ? 5 : 0;
    if (cfg->method == LC_SPEC_GB) params->gbEvalMode = GB_EVAL_GBLS;
    if (cfg->method == LC_SPEC_BLS) params->gbEvalMode = GB_EVAL_BLS;

    bool use_aov = (cfg->method == LC_SPEC_AOV);
    bool direct_grid = (cfg->method == LC_SPEC_GB || cfg->method == LC_SPEC_BLS);

    for (uint32_t i = 0; i < n; ++i) {
        buf->x[i] = data->x[i];
        buf->y[i] = data->y[i];
        buf->dy[i] = data->dy[i];
    }
    buf->n = n;
    preprocess_buffer(buf, params->epsilon, params->detrend_degree);

    double grid_span = buf->x[buf->n - 1] - buf->x[0];
    if (grid_span <= 0.0) grid_span = buf->x[buf->n - 1];
    const eval_method_t *eval_method = eval_method_for_params(params);
    double fmin = effective_grid_fmin(params, grid_span);
    double fstep;
    if (mode_uses_direct_eval_grid(params->mode)) {
        fstep = eval_method->direct_frequency_step(buf->n, grid_span, params);
    } else {
        fstep = 1.0 / ((double)params->nterms * (double)params->oversamplingFactor * grid_span * 0.5);
    }
    uint32_t nfreq = target_frequency_count(params, fmin, fstep);
    if (nfreq > buf->maxFreqCount) nfreq = buf->maxFreqCount;

    float *nll = (float *)calloc(nfreq > 0U ? (size_t)nfreq : 1U, sizeof(float));
    if (!nll) return -1;

    int result_code = -1;
    kt_forpool_t *pool = ctx->pool;

    if (direct_grid) {
        int rc = lc_run_direct_grid(buf, params, eval_method, pool, fmin, fstep, nfreq, nll, progress);
        if (rc == 1) { result_code = 1; goto done; }
        if (rc < 0) goto done;
    } else {
        if (!activate_target_nufft_plan(buf, params, nfreq)) goto done;

        bool ok;
        if (use_aov) {
            ok = lc_execute_aov_sweep_parallel(buf, params, pool, fmin, fstep, nfreq, progress, ctx->physical_cores);
        } else {
            ok = lc_execute_ihs_sweep_parallel(buf, params, pool, fmin, fstep, nfreq, progress, ctx->physical_cores);
        }
        if (!ok) goto done;
        if (lc_progress_is_cancelled(progress)) { result_code = 1; goto done; }

        int aov_n_eff = periodogram_effective_n(buf);
        int conv_threads = (pool && pool->n_threads > 1) ? pool->n_threads : 1;
        uint32_t conv_blocks = (uint32_t)conv_threads;
        if (conv_blocks > nfreq) conv_blocks = nfreq;
        lc_progress_add_total(progress, conv_blocks);
        lc_parallel_nll_convert(buf->power, nll, nfreq, params->nterms, aov_n_eff, use_aov, (params->statistic == STATISTIC_BAYES), pool, progress);
    }

    out->nfreq = nfreq;
    out->fmin = fmin;
    out->fstep = fstep;
    out->nll = nll;
    nll = NULL;

    lc_detect_peaks(out->nll, nfreq, fmin, fstep, cfg->oversampling > 0.0 ? cfg->oversampling : 5.0, cfg->peak_threshold, out->peaks, &out->npeaks);

    result_code = 0;

done:
    free(nll);
    return result_code;
}

/* Convenience wrapper: creates a context, runs one computation, destroys context. */
LC_API int lc_compute_periodogram(const lc_data_t *data, const lc_periodogram_config_t *cfg, lc_periodogram_result_t *out, lc_progress_t *progress) {
    int nthreads = (cfg && cfg->nthreads > 0) ? cfg->nthreads : 0;
    lc_compute_ctx_t *ctx = lc_compute_ctx_create(nthreads);
    if (!ctx) return -1;
    int rc = lc_compute_periodogram_ctx(ctx, data, cfg, out, progress);
    lc_compute_ctx_destroy(ctx);
    return rc;
}

/* Phased model overlay — included here (single-TU) to avoid duplicate symbols from
 * upstream header-only libraries (sds.h, convolution.h) that have non-static defs. */
#include "lc_model.c"
