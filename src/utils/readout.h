#ifndef READOUT_H
#define READOUT_H

#include <errno.h>
#include <fast_convert.h>
#include <fcntl.h>
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
    if (buffer->peaks) {
        free(buffer->peaks);
        buffer->peaks = NULL;
    }

    sdsfree(buffer->outBuf);
    sdsfree(buffer->spectrum);

    if (buffer->buf) {
        for (int i = 0; i < 3; i++) {
            if (buffer->buf && buffer->buf[i]) {
                free(buffer->buf[i]);
                buffer->buf[i] = NULL;
            }
        }
        if (buffer->buf) {
            free(buffer->buf);
            buffer->buf = NULL;
        }
    }
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
    bool needs_power_grid = !periodogram_uses_aov(params->periodogramMethod) || params->spectrum;
    if (needs_power_grid && !buffer->power) {
        buffer->power = aligned_alloc(64, round_buffer(((size_t)params->maxFreqCount + 2U) * sizeof(float)));
    }
    if (needs_power_grid && !buffer->power) goto error;
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
        buffer->buf = calloc(3, sizeof(void*));
        for (int i = 0; i < 3; i++) {
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

static inline void linreg_buffer(buffer_t* buffer) {
    // Center the measurement times - increases precision of future computation
    for (unsigned int i = 1; i < buffer->n; i++) {
        buffer->x[i] -= buffer->x[0];
    }
    buffer->x[0] = 0;

    // Initialize sums for regression
    double sumx = 0, sumxsq = 0, sumy = 0, sumxy = 0;
    double lin = 0, c = 0;

    for (unsigned int i = 0; i < buffer->n; i++) {
        sumx += buffer->x[i];                   // sum of x
        sumxsq += buffer->x[i] * buffer->x[i];  // sum of x^2
        sumy += buffer->y[i];                   // sum of y
        sumxy += buffer->x[i] * buffer->y[i];   // sum of x*y
    }

    // Calculate the denominator for the slope and intercept
    double denum = ((double)(buffer->n) * sumxsq) - (sumx * sumx);

    // Compute the slope (lin) and intercept (c)
    lin = ((buffer->n * sumxy) - (sumx * sumy)) / denum;
    c = ((sumy * sumxsq) - (sumx * sumxy)) / denum;

    // Adjust the y values based on the regression line
    for (unsigned int i = 0; i < buffer->n; i++) {
        buffer->y[i] -= (lin * buffer->x[i]) + c;
    }
    // sumy = 0;
    // for (unsigned int i = 0; i < buffer->n; i++) {sumy += buffer->y[i]};
    // sumy /= buffer->n;
    // for (unsigned int i = 0; i < buffer->n; i++) {buffer->y[i] -= (lin * buffer->x[i]) + c;}
}

static inline void linregw_buffer(buffer_t* buffer) {
    // Center the measurement times - increases precision of future computation
    for (unsigned int i = 1; i < buffer->n; i++) {
        buffer->x[i] -= buffer->x[0];
    }
    buffer->x[0] = 0;

    // Initialize sums for weighted regression
    double sumx = 0, sumxsq = 0, sumy = 0, sumxy = 0, sumw = 0;  // Add sum of weights
    double lin = 0, c = 0, w;

    for (unsigned int i = 0; i < buffer->n; i++) {
        w = 1 / (buffer->dy[i] * buffer->dy[i]);  // weight as the inverse of variance
        if (w > 0) {
            sumw += w;                                   // accumulate total weights
            sumx += buffer->x[i] * w;                    // weighted sum of x
            sumxsq += buffer->x[i] * buffer->x[i] * w;   // weighted sum of x^2
            sumy += buffer->y[i] * w;                    // weighted sum of y
            sumxy += (buffer->x[i] * buffer->y[i]) * w;  // weighted sum of x*y
        } else {
            buffer->dy[i] = 999.9;
        }
    }

    // Calculate the denominator for the slope and intercept
    double denum = (sumw * sumxsq) - (sumx * sumx);

    // Compute the slope (lin) and intercept (c)
    lin = ((sumw * sumxy) - (sumx * sumy)) / denum;
    c = ((sumy * sumxsq) - (sumx * sumxy)) / denum;
    buffer->magnitude = c;

    // Adjust the y values based on the regression line
    for (unsigned int i = 0; i < buffer->n; i++) {
        buffer->y[i] -= (lin * buffer->x[i]) + c;
    }
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

static inline void preprocess_buffer(buffer_t* buffer, double epsilon, int mode) {
    linregw_buffer(buffer);
    if (mode > 0) {
        linreg_buffer(buffer);
    }

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

    // Advise the kernel about sequential access
    posix_fadvise(fd, 0, file_size, POSIX_FADV_SEQUENTIAL);

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
