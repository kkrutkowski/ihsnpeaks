#ifndef READOUT_H
#define READOUT_H

#include <errno.h>
#include <fast_convert.h>
#include <fcntl.h>
#include <math.h>
#include <qfits/qfits.h>
#include <sds.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"

static inline size_t round_buffer(size_t size) { return (size + 63) & ~63; }

#include "convolution.h"

static inline void free_buffer(buffer_t* buffer) {
    buffer->allocated = false;
    if (buffer->x) {
        free(buffer->x);
        buffer->x = NULL;
    }
    if (buffer->y) {
        free(buffer->y);
        buffer->y = NULL;
    }
    if (buffer->dy) {
        free(buffer->dy);
        buffer->dy = NULL;
    }
    if (buffer->wy) {
        free(buffer->wy);
        buffer->wy = NULL;
    }
    if (buffer->pidx) {
        free(buffer->pidx);
        buffer->pidx = NULL;
    }
    if (buffer->readBuf) {
        free(buffer->readBuf);
        buffer->readBuf = NULL;
    }
    nufft1_free_workspace(buffer->nufftWorkspace);
    buffer->nufftWorkspace = NULL;
    free(buffer->power);
    buffer->power = NULL;
    buffer->powerCap = 0;
    free(buffer->blockReal);
    buffer->blockReal = NULL;
    free(buffer->blockImag);
    buffer->blockImag = NULL;
    free(buffer->inputReal);
    buffer->inputReal = NULL;
    free(buffer->inputImag);
    buffer->inputImag = NULL;
    free(buffer->workReal);
    buffer->workReal = NULL;
    free(buffer->workImag);
    buffer->workImag = NULL;
    free(buffer->deltaReal);
    buffer->deltaReal = NULL;
    free(buffer->deltaImag);
    buffer->deltaImag = NULL;
    free(buffer->deltaBackReal);
    buffer->deltaBackReal = NULL;
    free(buffer->deltaBackImag);
    buffer->deltaBackImag = NULL;
    free(buffer->fftReal);
    buffer->fftReal = NULL;
    free(buffer->fftImag);
    buffer->fftImag = NULL;
    free(buffer->cobraReal);
    buffer->cobraReal = NULL;
    free(buffer->cobraImag);
    buffer->cobraImag = NULL;
    free(buffer->aovScratch);
    buffer->aovScratch = NULL;
    buffer->aovScratchLen = 0;
    free(buffer->aovSw);
    buffer->aovSw = NULL;
    free(buffer->aovCw);
    buffer->aovCw = NULL;
    free(buffer->aovSyw);
    buffer->aovSyw = NULL;
    free(buffer->aovCyw);
    buffer->aovCyw = NULL;
    free(buffer->aovPower);
    buffer->aovPower = NULL;
    free(buffer->aovCondition);
    buffer->aovCondition = NULL;
    buffer->aovArrayCap = 0;
    buffer->aovTerms = 0;
    if (buffer->peaks) {
        free(buffer->peaks);
        buffer->peaks = NULL;
    }

    sdsfree(buffer->outBuf);
    sdsfree(buffer->spectrum);

    if (buffer->buf) {
        for (int i = 0; i < GB_SCRATCH_BUF_COUNT; i++) {
            if (buffer->buf[i]) {
                free(buffer->buf[i]);
                buffer->buf[i] = NULL;
            }
        }
        free(buffer->buf);
        buffer->buf = NULL;
    }
}

// Lazily allocate (or grow) AoV working arrays to exactly the required capacity.
static inline bool buffer_ensure_aov_arrays(buffer_t* buffer, size_t cap, int nterms) {
    if (buffer->aovSw && buffer->aovArrayCap >= cap && buffer->aovTerms >= nterms) return true;
    // Free old (too-small) arrays
    free(buffer->aovSw);
    buffer->aovSw = NULL;
    free(buffer->aovCw);
    buffer->aovCw = NULL;
    free(buffer->aovSyw);
    buffer->aovSyw = NULL;
    free(buffer->aovCyw);
    buffer->aovCyw = NULL;
    free(buffer->aovPower);
    buffer->aovPower = NULL;
    free(buffer->aovCondition);
    buffer->aovCondition = NULL;
    buffer->aovArrayCap = 0;
    buffer->aovTerms = 0;
    int max_factor = 2 * nterms;
    size_t trig_len = (size_t)(max_factor + 1) * cap;
    size_t ytrig_len = (size_t)(nterms + 1) * cap;
    buffer->aovSw = aligned_alloc(64, round_buffer(trig_len * sizeof(float)));
    buffer->aovCw = aligned_alloc(64, round_buffer(trig_len * sizeof(float)));
    buffer->aovSyw = aligned_alloc(64, round_buffer(ytrig_len * sizeof(float)));
    buffer->aovCyw = aligned_alloc(64, round_buffer(ytrig_len * sizeof(float)));
    buffer->aovPower = aligned_alloc(64, round_buffer(cap * sizeof(float)));
    buffer->aovCondition = aligned_alloc(64, round_buffer(cap * sizeof(float)));
    if (!buffer->aovSw || !buffer->aovCw || !buffer->aovSyw || !buffer->aovCyw || !buffer->aovPower || !buffer->aovCondition) return false;
    buffer->aovArrayCap = cap;
    buffer->aovTerms = nterms;
    return true;
}

// Lazily allocate (or grow) the power grid to exactly the required capacity.
// Sweep-time callers pass slice_nfreq+1 so parallel workers only allocate what
// their frequency slice actually needs; the primary gets the full grid up front.
static inline bool buffer_ensure_power(buffer_t* buffer, size_t cap) {
    if (buffer->power && buffer->powerCap >= cap) return true;
    float* np = aligned_alloc(64, round_buffer(cap * sizeof(float)));
    if (!np) return false;
    free(buffer->power);
    buffer->power = np;
    buffer->powerCap = cap;
    return true;
}

static inline int alloc_buffer(buffer_t* buffer, parameters* params) {
    buffer->len = params->maxLen;
    buffer->allocated = true;
    buffer->terms = params->nterms;
    buffer->maxFreqCount = params->maxFreqCount;
    buffer->paddedLen = (uint32_t)(((size_t)params->maxLen + (size_t)VEC_LEN - 1U) & ~((size_t)VEC_LEN - 1U));
    uint32_t max_plan_idx = params->nufftPlanCount > 0U ? params->nufftPlanCount - 1U : 0U;
    if (params->nufftPlanCount > 0U) {
        buffer->activePlanIndex = max_plan_idx;
        buffer->activeGridLen = params->nufftPlanCache[max_plan_idx].gridLen;
        buffer->activeOutputLen = params->nufftPlanCache[max_plan_idx].outputLen;
        buffer->activeLadderLevels = params->ladderLevels;
    }
    buffer->spectrum = sdsempty();
    buffer->outBuf = sdsempty();
    if (!buffer->x) {
        buffer->x = aligned_alloc(64, round_buffer(params->maxLen * sizeof(double)));
    }
    if (!buffer->x) goto error;
    if (!buffer->y) {
        buffer->y = aligned_alloc(64, round_buffer(params->maxLen * sizeof(float)));
    }
    if (!buffer->y) goto error;
    if (!buffer->dy) {
        buffer->dy = aligned_alloc(64, round_buffer(params->maxLen * sizeof(float)));
    }
    if (!buffer->dy) goto error;
    if (!buffer->wy) {
        buffer->wy = aligned_alloc(64, round_buffer(params->maxLen * sizeof(float)));
    }
    if (!buffer->wy) goto error;
    if (!buffer->pidx) {
        buffer->pidx = (size_t*)malloc(1024 * sizeof(size_t));
    }
    if (!buffer->pidx) goto error;
    if (!buffer->readBuf) {
        buffer->readBuf = aligned_alloc(64, round_buffer(params->maxSize));
    }
    if (!buffer->readBuf) goto error;
    if (!buffer->power) {
        // Power grid: full frequency grid only when spectrum mode needs it;
        // otherwise lazily grown per slice by buffer_ensure_power() at sweep time.
        bool needs_power_grid = mode_uses_direct_eval_grid(params->mode) || !periodogram_uses_aov(params->periodogramMethod) || params->spectrum;
        if (needs_power_grid) {
            size_t power_entries = params->spectrum ? ((size_t)params->maxFreqCount + 2U) : ((size_t)params->outputLen + 2U);
            buffer->power = aligned_alloc(64, round_buffer(power_entries * sizeof(float)));
            if (buffer->power) buffer->powerCap = power_entries;
        }
    }
    if (params->spectrum && !buffer->power) goto error;
    if (mode_uses_direct_eval_grid(params->mode) && !buffer->power) goto error;
    if (!periodogram_uses_aov(params->periodogramMethod) && !buffer->power) goto error;
    if (!buffer->blockReal) {
        buffer->blockReal = aligned_alloc(64, round_buffer((size_t)params->outputLen * sizeof(float)));
    }
    if (!buffer->blockReal) goto error;
    if (!buffer->blockImag) {
        buffer->blockImag = aligned_alloc(64, round_buffer((size_t)params->outputLen * sizeof(float)));
    }
    if (!buffer->blockImag) goto error;
    if (!buffer->inputReal) {
        buffer->inputReal = aligned_alloc(64, round_buffer((size_t)buffer->paddedLen * sizeof(float)));
    }
    if (!buffer->inputReal) goto error;
    if (!buffer->inputImag) {
        buffer->inputImag = aligned_alloc(64, round_buffer((size_t)buffer->paddedLen * sizeof(float)));
    }
    if (!buffer->inputImag) goto error;
    size_t ladder_len = (size_t)NUFFT_LADDER_LEVEL_CAP * (size_t)buffer->paddedLen;
    if (!buffer->workReal) {
        buffer->workReal = aligned_alloc(64, round_buffer(ladder_len * sizeof(float)));
    }
    if (!buffer->workReal) goto error;
    if (!buffer->workImag) {
        buffer->workImag = aligned_alloc(64, round_buffer(ladder_len * sizeof(float)));
    }
    if (!buffer->workImag) goto error;
    if (!buffer->deltaReal) {
        buffer->deltaReal = aligned_alloc(64, round_buffer(ladder_len * sizeof(float)));
    }
    if (!buffer->deltaReal) goto error;
    if (!buffer->deltaImag) {
        buffer->deltaImag = aligned_alloc(64, round_buffer(ladder_len * sizeof(float)));
    }
    if (!buffer->deltaImag) goto error;
    if (!buffer->deltaBackReal) {
        buffer->deltaBackReal = aligned_alloc(64, round_buffer(ladder_len * sizeof(float)));
    }
    if (!buffer->deltaBackReal) goto error;
    if (!buffer->deltaBackImag) {
        buffer->deltaBackImag = aligned_alloc(64, round_buffer(ladder_len * sizeof(float)));
    }
    if (!buffer->deltaBackImag) goto error;
    if (!buffer->fftReal) {
        buffer->fftReal = aligned_alloc(64, round_buffer(params->nufftExternalSizes.fft_len * sizeof(float)));
    }
    if (!buffer->fftReal) goto error;
    if (!buffer->fftImag) {
        buffer->fftImag = aligned_alloc(64, round_buffer(params->nufftExternalSizes.fft_len * sizeof(float)));
    }
    if (!buffer->fftImag) goto error;
    if (!buffer->cobraReal) {
        buffer->cobraReal = aligned_alloc(64, round_buffer(params->nufftExternalSizes.cobra_len * sizeof(float)));
    }
    if (!buffer->cobraReal) goto error;
    if (!buffer->cobraImag) {
        buffer->cobraImag = aligned_alloc(64, round_buffer(params->nufftExternalSizes.cobra_len * sizeof(float)));
    }
    if (!buffer->cobraImag) goto error;
    if (!buffer->nufftWorkspace) {
        if (params->nufftPlanCount == 0U) goto error;
        buffer->nufftWorkspace =
            nufft1_create_workspace(params->nufftPlanCache[max_plan_idx].plan, buffer->fftReal, buffer->fftImag, buffer->cobraReal, buffer->cobraImag);
    }
    if (!buffer->nufftWorkspace) goto error;
    if (periodogram_uses_aov(params->periodogramMethod) && !buffer->aovScratch) {
        buffer->aovScratchLen = params->nterms > 0 ? (size_t)(26 * params->nterms + 16) : 0U;
        if (buffer->aovScratchLen > 0U) {
            buffer->aovScratch = aligned_alloc(64, round_buffer(buffer->aovScratchLen * sizeof(float)));
        }
        if (!buffer->aovScratch) goto error;
    }
    if (!buffer->peaks) {
        size_t peak_count = params->prewhiten ? (size_t)2 * (size_t)params->npeaks : (size_t)params->npeaks;
        buffer->peaks = calloc(peak_count, sizeof(peak_t));
    }
    buffer->nPeaks = 0;
    if (!buffer->peaks) goto error;
    if ((params->mode > 0 || params->prewhiten) && !buffer->buf) {
        buffer->buf = calloc(GB_SCRATCH_BUF_COUNT, sizeof(void*));
        for (int i = 0; i < GB_SCRATCH_BUF_COUNT; i++) {
            buffer->buf[i] = aligned_alloc(64, round_buffer(params->maxLen * sizeof(uint64_t)));
        }
    }
    return 0;

error:
    free_buffer(buffer);
    fprintf(stderr, "Failed to allocate buffer\n");
    return -1;
}

static inline void print_buffer(buffer_t* buffer) {
    for (int i = 0; i < buffer->n; i++) {
        printf("%.2f\t%.2f\t%.2f\n", buffer->x[i], buffer->y[i], buffer->dy[i]);
    }
}

static inline void detrend_buffer_szego(buffer_t* buffer, double epsilon, int detrend_degree) {
    if (!buffer || buffer->n == 0) return;
    uint32_t n = buffer->n;

    // Center the measurement times - increases precision of future computation
    double x0 = buffer->x[0];
    for (uint32_t i = 1; i < n; i++) {
        buffer->x[i] -= x0;
    }
    buffer->x[0] = 0.0;

    double span = buffer->x[n - 1] - buffer->x[0];

    // Compute weights and weighted mean
    double ws = 0.0;
    double wysum = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        double w = 1.0;
        if (buffer->dy) {
            double dy_val = (double)buffer->dy[i];
            if (dy_val > 0.0 && isfinite(dy_val)) {
                w = 1.0 / (dy_val * dy_val + epsilon);
            } else {
                buffer->dy[i] = 999.9f;
                w = 1.0 / (999.9 * 999.9 + epsilon);
            }
        }
        ws += w;
        wysum += w * (double)buffer->y[i];
    }

    double y_mean = (ws > 0.0 && isfinite(ws)) ? (wysum / ws) : 0.0;
    buffer->magnitude = (float)y_mean;

    // Degrees of freedom clamping: 4d + 2 <= n -> d <= (n - 2) / 4
    int req_degree = detrend_degree < 0 ? 0 : detrend_degree;
    int max_supported_degree = ((int)n - 2) / 4;
    if (max_supported_degree < 0) max_supported_degree = 0;
    int eff_degree = req_degree < max_supported_degree ? req_degree : max_supported_degree;

    if (!(span > 0.0) || !isfinite(span) || n < 6 || eff_degree == 0) {
        // Degree 0: Subtract weighted mean only
        for (uint32_t i = 0; i < n; i++) {
            buffer->y[i] -= (float)y_mean;
        }
        return;
    }

    int nh = eff_degree;  // number of half-integer harmonics at f_base = 1 / (2 * span)
    int nn2 = 2 * nh;

    double* rw = (double*)malloc(n * sizeof(double));
    double* zr = (double*)malloc(n * sizeof(double));
    double* zi = (double*)malloc(n * sizeof(double));
    double* pr = (double*)malloc(n * sizeof(double));
    double* pi = (double*)malloc(n * sizeof(double));
    double* znr = (double*)malloc(n * sizeof(double));
    double* zni = (double*)malloc(n * sizeof(double));
    double* c_nh = (double*)malloc(n * sizeof(double));
    double* s_nh = (double*)malloc(n * sizeof(double));
    double* cfr = (double*)malloc(n * sizeof(double));
    double* cfi = (double*)malloc(n * sizeof(double));
    double* my = (double*)calloc(n, sizeof(double));

    if (!rw || !zr || !zi || !pr || !pi || !znr || !zni || !c_nh || !s_nh || !cfr || !cfi || !my) {
        free(rw);
        free(zr);
        free(zi);
        free(pr);
        free(pi);
        free(znr);
        free(zni);
        free(c_nh);
        free(s_nh);
        free(cfr);
        free(cfi);
        free(my);
        for (uint32_t i = 0; i < n; i++) buffer->y[i] -= (float)y_mean;
        return;
    }

    double f_base = 0.5 / span;
    double t_mid = span / 2.0;

    for (uint32_t i = 0; i < n; i++) {
        double w = 1.0;
        if (buffer->dy) {
            double dy_val = (double)buffer->dy[i];
            if (dy_val > 0.0 && isfinite(dy_val)) {
                w = 1.0 / (dy_val * dy_val + epsilon);
            } else {
                w = 1.0 / (999.9 * 999.9 + epsilon);
            }
        }
        rw[i] = sqrt(w);

        double phase = (buffer->x[i] - t_mid) * f_base;
        zr[i] = cos(2.0 * M_PI * phase);
        zi[i] = sin(2.0 * M_PI * phase);
        pr[i] = rw[i];
        pi[i] = 0.0;
        znr[i] = 1.0;
        zni[i] = 0.0;

        double angle_nh = 2.0 * M_PI * (double)nh * phase;
        double cnh = cos(angle_nh);
        double snh = sin(angle_nh);
        c_nh[i] = cnh;
        s_nh[i] = snh;

        double vy = ((double)buffer->y[i] - y_mean) * rw[i];
        cfr[i] = vy * cnh;
        cfi[i] = vy * snh;
    }

    for (int step = 0; step <= nn2; step++) {
        double sn = 0.0;
        double alr = 0.0, ali = 0.0;
        double scr = 0.0, sci = 0.0;

        for (uint32_t i = 0; i < n; i++) {
            double pr_i = pr[i];
            double pi_i = pi[i];
            sn += pr_i * pr_i + pi_i * pi_i;
            alr += (zr[i] * pr_i - zi[i] * pi_i) * rw[i];
            ali += (zr[i] * pi_i + zi[i] * pr_i) * rw[i];

            scr += pr_i * cfr[i] + pi_i * cfi[i];
            sci += pr_i * cfi[i] - pi_i * cfr[i];
        }

        if (sn < 1e-30) sn = 1e-30;
        alr /= sn;
        ali /= sn;

        double cr = scr / sn;
        double ci = sci / sn;

        for (uint32_t i = 0; i < n; i++) {
            double pr_i = pr[i];
            double pi_i = pi[i];
            double Ay = cr * pr_i - ci * pi_i;
            double By = cr * pi_i + ci * pr_i;
            my[i] += Ay * c_nh[i] + By * s_nh[i];
        }

        for (uint32_t i = 0; i < n; i++) {
            double pr_i = pr[i];
            double pi_i = pi[i];
            double zr_i = zr[i];
            double zi_i = zi[i];
            double znr_i = znr[i];
            double zni_i = zni[i];

            double sr = alr * znr_i - ali * zni_i;
            double si = alr * zni_i + ali * znr_i;

            double new_pr = pr_i * zr_i - pi_i * zi_i - sr * pr_i - si * pi_i;
            double new_pi = pr_i * zi_i + pi_i * zr_i + sr * pi_i - si * pr_i;
            pr[i] = new_pr;
            pi[i] = new_pi;

            double new_znr = znr_i * zr_i - zni_i * zi_i;
            double new_zni = zni_i * zr_i + znr_i * zi_i;
            znr[i] = new_znr;
            zni[i] = new_zni;
        }
    }

    for (uint32_t i = 0; i < n; i++) {
        double P_y = y_mean + (rw[i] > 0.0 ? my[i] / rw[i] : 0.0);
        buffer->y[i] = (float)((double)buffer->y[i] - P_y);
    }

    free(rw);
    free(zr);
    free(zi);
    free(pr);
    free(pi);
    free(znr);
    free(zni);
    free(c_nh);
    free(s_nh);
    free(cfr);
    free(cfi);
    free(my);
}

static inline void refresh_weighted_signal_buffer(buffer_t* buffer, double epsilon) {
    double wysum = 0.0;
    double wysqsum = 0.0;
    for (uint32_t i = 0; i < buffer->n; i++) {
        float weight = 1.0f / ((buffer->dy[i] * buffer->dy[i]) + (float)epsilon);
        buffer->wy[i] = weight * buffer->y[i];
        wysum += fabs(buffer->wy[i]);
        wysqsum += buffer->wy[i] * buffer->wy[i];
    }
    buffer->amp_neff = wysqsum > 0.0 ? ((wysum * wysum) / wysqsum) - 2.0 : 0.0;
    double norm_neff = buffer->amp_neff > 0.0 ? buffer->amp_neff : buffer->neff;
    double norm = (wysum > 0.0 && norm_neff > 0.0 && double_is_finite_bits(wysum) && double_is_finite_bits(norm_neff)) ? sqrt(norm_neff) / wysum : 0.0;
    if (!double_is_finite_bits(norm)) norm = 0.0;
    for (uint32_t i = 0; i < buffer->n; i++) {
        buffer->wy[i] *= (float)norm;
    }
}

static inline void preprocess_buffer(buffer_t* buffer, double epsilon, int detrend_degree) {
    detrend_buffer_szego(buffer, epsilon, detrend_degree);

    double wsum = 0;
    double wsqsum = 0;

    for (uint32_t i = 0; i < buffer->n; i++) {
        float weight = 1.0f / ((buffer->dy[i] * buffer->dy[i]) + (float)epsilon);
        wsum += fabs(weight);
        wsqsum += weight * weight;
    }
    buffer->neff = wsqsum > 0.0 ? ((wsum * wsum) / wsqsum) - 2.0 : 0.0;
    refresh_weighted_signal_buffer(buffer, epsilon);
}

static inline void read_fits(const char* in_file, buffer_t* buffer) {
    qfits_table* table = qfits_table_open(in_file, 1);
    if (!table) {
        fprintf(stderr, "Failed to open FITS table in %s\n", in_file);
        buffer->n = 0;
        return;
    }

    int col_time = -1, col_flux = -1, col_fluxerr = -1;
    for (int i = 0; i < table->nc; i++) {
        if (strcmp(table->col[i].tlabel, "TIME") == 0)
            col_time = i;
        else if (strstr(table->col[i].tlabel, "PDCSAP_FLUX") != NULL || strcmp(table->col[i].tlabel, "SAP_FLUX") == 0 ||
                 strcmp(table->col[i].tlabel, "FLUX") == 0)
            col_flux = i;
        else if (strstr(table->col[i].tlabel, "PDCSAP_FLUX_ERR") != NULL || strcmp(table->col[i].tlabel, "SAP_FLUX_ERR") == 0 ||
                 strcmp(table->col[i].tlabel, "FLUX_ERR") == 0)
            col_fluxerr = i;
    }

    if (col_time < 0 || col_flux < 0 || col_fluxerr < 0) {
        fprintf(stderr, "Failed to find required columns (TIME, FLUX, FLUX_ERR) in %s\n", in_file);
        qfits_table_close(table);
        buffer->n = 0;
        return;
    }

    double null_time = 0.0, null_flux = 0.0, null_fluxerr = 0.0;
    void* time_data = qfits_query_column_data(table, col_time, NULL, &null_time);
    void* flux_data = qfits_query_column_data(table, col_flux, NULL, &null_flux);
    void* fluxerr_data = qfits_query_column_data(table, col_fluxerr, NULL, &null_fluxerr);

    if (!time_data || !flux_data || !fluxerr_data) {
        fprintf(stderr, "Failed to read column data from %s\n", in_file);
        free(time_data);
        free(flux_data);
        free(fluxerr_data);
        qfits_table_close(table);
        buffer->n = 0;
        return;
    }

    tfits_type t_type = table->col[col_time].atom_type;
    tfits_type f_type = table->col[col_flux].atom_type;
    tfits_type fe_type = table->col[col_fluxerr].atom_type;

    size_t idx = 0;
    for (int i = 0; i < table->nr && idx < buffer->len; i++) {
        double t = 0.0, f = 0.0, fe = 0.0;

        if (t_type == TFITS_BIN_TYPE_D)
            t = ((double*)time_data)[i];
        else if (t_type == TFITS_BIN_TYPE_E)
            t = ((float*)time_data)[i];

        if (f_type == TFITS_BIN_TYPE_D)
            f = ((double*)flux_data)[i];
        else if (f_type == TFITS_BIN_TYPE_E)
            f = ((float*)flux_data)[i];

        if (fe_type == TFITS_BIN_TYPE_D)
            fe = ((double*)fluxerr_data)[i];
        else if (fe_type == TFITS_BIN_TYPE_E)
            fe = ((float*)fluxerr_data)[i];

        if (isnan(t) || isinf(t) || isnan(f) || isinf(f) || isnan(fe) || isinf(fe)) continue;
        if (t == null_time || f == null_flux || fe == null_fluxerr) continue;

        buffer->x[idx] = t;
        buffer->y[idx] = (float)f;
        buffer->dy[idx] = (float)fe;
        idx++;
    }

    buffer->n = idx;

    free(time_data);
    free(flux_data);
    free(fluxerr_data);
    qfits_table_close(table);
}

static inline void read_dat(const char* in_file, buffer_t* buffer) {
    // Open the file
    int fd = open(in_file, O_RDONLY);
    if (fd == -1) {
        perror("Error opening file");
        // return;
    }

    // Get the size of the file
    struct stat file_stat;
    if (fstat(fd, &file_stat) == -1) {
        perror("Failed to get file size");
        close(fd);
        // return;
    }
    size_t file_size = file_stat.st_size;

    // Advise the kernel about sequential access when the platform exposes it.
#if defined(__linux__) && defined(POSIX_FADV_SEQUENTIAL)
    posix_fadvise(fd, 0, file_size, POSIX_FADV_SEQUENTIAL);
#endif

    // Use existing buffer to avoid reallocation
    char* dataBuffer = buffer->readBuf;

    // Read the file in one go
    ssize_t bytes_read = read(fd, dataBuffer, file_size);
    if (bytes_read == -1) {
        perror("Failed to read file");
        close(fd);
        return;
    } else if (bytes_read == 0) {
        fprintf(stderr, "No data read from file\n");
        close(fd);
        return;
    }

    // Null-terminate the buffer and close the file
    close(fd);
    dataBuffer[bytes_read] = '\0';

    // Parse the data
    char* it = dataBuffer;
    char* end = dataBuffer + bytes_read;
    double tempX;
    float tempY, tempDY;
    size_t idx = 0;

    while (it < end && idx < buffer->len) {
        // Parse temporary variables
        tempX = fast_strtod(it, &it);
        if (it == NULL || it >= end) break;
        it++;
        tempY = fast_strtof(it, &it);
        if (it == NULL || it >= end) break;
        it++;
        tempDY = fast_strtof(it, &it);
        if (it == NULL || it >= end) break;
        it++;

        // printf("%f\t%f\t%f\n", tempX, tempY, tempDY); // ok
        buffer->x[idx] = tempX;
        buffer->y[idx] = tempY;
        buffer->dy[idx] = tempDY;
        idx++;
    }

    buffer->n = idx;
}

static inline void append_peak(buffer_t* buff, const parameters* params, const double freq, const float magnitude, double df) {
    const eval_method_t* method = eval_method_for_params(params);
    peak_t appended = {0};  // peak_t tmp;
    appended.freq = freq;
    appended.p = magnitude;
    const int maxPeaks = params->npeaks;
    const int mode = params->mode;
    int idx = buff->nPeaks;
    if (mode_defers_peak_evaluation(mode)) {
        while (idx > 0 && magnitude > buff->peaks[idx - 1].p) {
            idx--;
        }
        if (buff->nPeaks < maxPeaks) {
            buff->nPeaks++;
        }
        for (int i = buff->nPeaks - 1; i > idx; i--) {
            buff->peaks[i] = buff->peaks[i - 1];
        }
        if (idx < maxPeaks) {
            buff->peaks[idx] = appended;
        }
    } else {
        if (mode_eagerly_refines_peaks(mode)) {
            binsearch_peak(&appended, buff, params, method, df);
        } else {
            evaluate_peak_at_current_frequency(&appended, buff, params, method);
        }
        float rank = eval_peak_rank(method, &appended);
        while (idx > 0 && rank > eval_peak_rank(method, &buff->peaks[idx - 1])) {
            idx--;
        }
        if (buff->nPeaks < maxPeaks) {
            buff->nPeaks++;
        }
        for (int i = buff->nPeaks - 1; i > idx; i--) {
            buff->peaks[i] = buff->peaks[i - 1];
        }
        if (idx < maxPeaks) {
            buff->peaks[idx] = appended;
        }
    }
}

#endif
