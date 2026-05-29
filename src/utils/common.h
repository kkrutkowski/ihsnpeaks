#ifndef COMMON_H
#define COMMON_H

#include <klib/kvec.h>
#include <pthread.h>
#include <sds.h>
#include <stdbool.h>
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

typedef enum { GB_EVAL_GBLS = 0, GB_EVAL_GBAS } gb_eval_mode;

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

static inline const char* gb_eval_name(gb_eval_mode mode) { return mode == GB_EVAL_GBAS ? "gbas" : "gbls"; }

static inline const char* gb_stat_label(gb_eval_mode mode) { return mode == GB_EVAL_GBAS ? "T" : "R2"; }

typedef struct {
    double freq;
    float p;
    float amp;
    float r2;
} peak_t;

typedef struct {
    bool allocated;
    uint8_t loc_iter;
    uint32_t len;
    uint32_t n;
    uint32_t r2;
    uint32_t terms;
    uint32_t memBlockSize;
    double neff;

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
    nufft1_workspace* nufftWorkspace;

    sds spectrum;
    sds outBuf;
} buffer_t;

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
    nufft1_plan* nufftPlan;
    float* nufftTwiddleReal;
    float* nufftTwiddleImag;

    kvec_target_t targets;

    buffer_t** buffers;
    int nbuffers;
} parameters;

#endif  // COMMON_H
