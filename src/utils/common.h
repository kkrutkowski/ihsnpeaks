#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <complex.h>

#include <sds.h>
#include <klib/kvec.h>
#include <fftw3.h>

typedef struct {
    double freq;
    float  p;
    float  amp;
    float  r;
} peak_t;

typedef struct {
    bool allocated;
    uint8_t loc_iter;
    uint32_t   len;
    uint32_t     n;
    uint32_t     r;
    uint32_t terms;
    uint32_t memBlockSize;

    peak_t* peaks;
    uint32_t nPeaks;

    char*  readBuf;
    double*      x;
    float*       y;
    float*      dy;

    size_t* pidx; //phase indices for counting sort

    void** buf;

    float magnitude;
    uint32_t gridSize; uint32_t nGrids;

    complex float ** grids; //used to compute the FFT
    uint32_t** gidx; //grid indices for NFFT
    float** gdist;
    float** weights;

    sds spectrum;
    sds outBuf;
} buffer_t;

typedef union {
    uint64_t data; // The raw 64-bit representation
    struct {
        uint16_t key;        // 16-bit key
        uint16_t idx; // 16-bit index
        float val;         // 32-bit floating-point value
    } parts;
} kvpair;


// Define the target struct to hold file path and
typedef struct {char *path;} target_t;
typedef kvec_t(target_t) kvec_target_t;


typedef struct {
    char*  target;
    char* outFile;
    float fmin;    float fmax;
    float threshold;
    float oversamplingFactor;
    float epsilon;
    int npeaks;    int nterms;
    int mode; int jobs;
    bool isFile;
    bool spectrum;
    bool debug;
    bool corrected;
    bool idle;
    bool prewhiten;
    pthread_mutex_t mutex; pthread_mutex_t counter_mutex; int iter_count;

    fftwf_plan plan;

    //Variables used to estimate optimal hyperparameters from metadata
    uint32_t gridRatio; uint32_t defaultGridRatio;
    //uint32_t blockSize;       //size of the FFT plan applied
    //uint32_t bufferSize;      //size of the FFT buffer
    uint32_t maxSize;         //number of bytes in the longest time series of the processed batch
    uint32_t maxLen;          //number of measurements in the longest time series of the processed batch
    uint32_t gridLen;         //length of the FFT grid used in transform
    uint64_t avgLen;          //average number of measurements per time series in the batch

    kvec_target_t targets;

    buffer_t** buffers; int nbuffers;
} parameters;

#endif //COMMON_H
