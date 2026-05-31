#ifndef COMMON_H
#define COMMON_H

#include <klib/kvec.h>
#include <math.h>
#include <pthread.h>
#include <sds.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "../nufft/nufft1.h"

typedef enum {
    PERIODOGRAM_IHS = 0,
    PERIODOGRAM_AOV,
    PERIODOGRAM_AOVMH,
    PERIODOGRAM_AOBMHW,
    PERIODOGRAM_CHI,
    PERIODOGRAM_CHI2,
    PERIODOGRAM_FASTCHI2
} periodogram_method;

typedef enum { GB_EVAL_GBLS = 0, GB_EVAL_GBAW } gb_eval_mode;

static inline bool float_is_nan_bits(float value) {
    union {
        float f;
        uint32_t u;
    } bits = {value};
    return (bits.u & UINT32_C(0x7fffffff)) > UINT32_C(0x7f800000);
}

static inline bool float_is_finite_bits(float value) {
    union {
        float f;
        uint32_t u;
    } bits = {value};
    return (bits.u & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

static inline bool double_is_finite_bits(double value) {
    union {
        double f;
        uint64_t u;
    } bits = {value};
    return (bits.u & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}

static inline bool periodogram_uses_aov(periodogram_method method) { return method != PERIODOGRAM_IHS; }

static inline const char* periodogram_method_name(periodogram_method method) {
    switch (method) {
        case PERIODOGRAM_AOV:
            return "aov";
        case PERIODOGRAM_AOVMH:
            return "aovmh";
        case PERIODOGRAM_AOBMHW:
            return "aobmhw";
        case PERIODOGRAM_CHI:
            return "chi";
        case PERIODOGRAM_CHI2:
            return "chi2";
        case PERIODOGRAM_FASTCHI2:
            return "fastchi2";
        case PERIODOGRAM_IHS:
        default:
            return "ihs";
    }
}

static inline const char* gb_eval_name(gb_eval_mode mode) { return mode == GB_EVAL_GBAW ? "gbaw" : "gbls"; }

static inline const char* gb_stat_label(gb_eval_mode mode) { return mode == GB_EVAL_GBAW ? "T" : "R2"; }

static inline bool mode_uses_direct_gb_grid(int mode) { return mode >= 5; }

static inline bool mode_defers_peak_evaluation(int mode) { return mode < 3 || mode == 5; }

static inline bool mode_refines_retained_peaks(int mode) { return mode == 2 || mode == 3 || mode == 5; }

static inline bool mode_eagerly_refines_peaks(int mode) { return mode == 4 || mode == 6; }

static inline bool mode_evaluates_all_local_peaks(int mode) { return mode == 3 || mode == 4 || mode == 6; }

static inline int gb_convolution_radius(uint32_t n, float gbAlpha) {
    int r = (int)ceil(0.5 * sqrt((2.0 * (double)n) + 1.0) + ((double)n * 0.25 * (double)gbAlpha));
    r >>= 1;
    if (r < 1) r = 1;
    return r;
}

static inline double gb_direct_frequency_step(uint32_t n, double time_span, float oversamplingFactor, float gbAlpha) {
    if (n == 0 || time_span <= 0.0 || oversamplingFactor <= 0.0f) return 0.0;
    int r = gb_convolution_radius(n, gbAlpha);
    return (12.0 * (double)r) / ((double)n * (double)oversamplingFactor * time_span);
}

typedef struct {
    double freq;
    float p;
    float amp;
    float r2;
} peak_t;

enum { NUFFT_LADDER_LEVEL_CAP = 8U };

typedef struct {
    uint32_t gridLen;
    uint32_t outputLen;
    nufft1_external_sizes externalSizes;
    nufft1_plan* plan;
    float* twiddleReal;
    float* twiddleImag;
} nufft_plan_cache_entry_t;

typedef struct {
    bool allocated;
    uint8_t loc_iter;
    uint32_t len;
    uint32_t n;
    uint32_t r2;
    uint32_t terms;
    uint32_t memBlockSize;
    double neff;
    double amp_neff;

    peak_t* peaks;
    uint32_t nPeaks;

    char* readBuf;
    double* x;
    float* y;
    float* dy;
    float* wy;

    size_t* pidx;  // phase indices for counting sort

    void** buf;

    float magnitude;
    uint32_t maxFreqCount;
    uint32_t paddedLen;

    float* power;
    float* blockReal;
    float* blockImag;
    float* inputReal;
    float* inputImag;
    float* workReal;
    float* workImag;
    float* deltaReal;
    float* deltaImag;
    float* fftReal;
    float* fftImag;
    float* cobraReal;
    float* cobraImag;
    float* aovScratch;
    size_t aovScratchLen;
    nufft1_workspace* nufftWorkspace;
    uint32_t activePlanIndex;
    uint32_t activeGridLen;
    uint32_t activeOutputLen;
    uint32_t activeLadderLevels;

    sds spectrum;
    sds outBuf;
} buffer_t;

static inline int periodogram_effective_n(const buffer_t* buffer) {
    if (!buffer || buffer->neff <= 0.0) return 0;
    double n_eff = floor(buffer->neff);
    if (n_eff < 1.0) return 0;
    if (n_eff > (double)INT32_MAX) return INT32_MAX;
    return (int)n_eff;
}

typedef union {
    uint64_t data;  // The raw 64-bit representation
    struct {
        uint16_t key;  // 16-bit key
        uint16_t idx;  // 16-bit index
        float val;     // 32-bit floating-point value
    } parts;
} kvpair;

// Define the target struct to hold file path and
typedef struct {
    char* path;
} target_t;
typedef kvec_t(target_t) kvec_target_t;

typedef struct {
    char* target;
    char* outFile;
    float fmin;
    float fmax;
    float threshold;
    float r2_threshold;
    float oversamplingFactor;
    float epsilon;
    float gbAlpha;
    int npeaks;
    int nterms;
    int mode;
    int jobs;
    periodogram_method periodogramMethod;
    gb_eval_mode gbEvalMode;
    bool isFile;
    bool spectrum;
    bool debug;
    bool corrected;
    bool idle;
    bool prewhiten;
    bool outputPeriod;
    pthread_mutex_t mutex;
    pthread_mutex_t counter_mutex;
    int iter_count;

    // Variables used to estimate optimal hyperparameters from metadata
    uint32_t gridRatio;
    uint32_t defaultGridRatio;
    uint32_t maxSize;  // number of bytes in the longest time series of the processed batch
    uint32_t maxLen;   // number of measurements in the longest time series of the processed batch
    uint32_t gridLen;  // length of the NuFFT block plan
    uint32_t outputLen;
    uint32_t maxFreqCount;
    uint32_t ladderLevels;
    uint64_t avgLen;  // average number of measurements per time series in the batch
    double maxTimeSpan;
    nufft1_mode gridMode;
    nufft1_external_sizes nufftExternalSizes;
    nufft_plan_cache_entry_t* nufftPlanCache;
    float* nufftTwiddleReal;
    float* nufftTwiddleImag;
    uint32_t nufftPlanCount;

    kvec_target_t targets;

    buffer_t** buffers;
    int nbuffers;
} parameters;

#endif  // COMMON_H
