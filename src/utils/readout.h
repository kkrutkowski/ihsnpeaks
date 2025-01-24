#ifndef READOUT_H
#define READOUT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <complex.h>

#include <sds.h>
#include <fast_convert.h>
#include <fftw3.h>

//#include "../include/mufft/mufft.x86.h"

static inline size_t round_buffer(size_t size) {return (size + 63) & ~63;}

typedef struct {
    double freq;
    float  p;
    float  amp;
    float  chi2;
} peak_t;

typedef struct {
    bool allocated;
    uint8_t loc_iter;
    uint32_t   len;
    uint32_t     n;
    uint32_t terms;
    uint32_t memBlockSize;

    peak_t* peaks;
    uint32_t nPeaks;

    char*  readBuf;
    double*      x;
    float*       y;
    float*      dy;

    uint16_t* pidx; //phase indices for counting sort

    //D_VEC*    bufd1;
    //D_VEC*    bufd2;

    float magnitude;
    uint32_t gridSize; uint32_t nGrids;

    complex float ** grids; //used to compute the FFT
    uint32_t** gidx; //grid indices for NFFT
    float** gdist;
    float** weights;

    sds spectrum;
    sds outBuf;
} buffer_t;

static inline void free_buffer (buffer_t* buffer) {
    buffer -> allocated = false;
    if (buffer -> x)       {free(buffer -> x);  buffer -> x = NULL;}
    if (buffer -> y)       {free(buffer -> y);  buffer -> y = NULL;}
    if (buffer -> dy)      {free(buffer -> dy); buffer-> dy  = NULL;}
    if (buffer -> pidx)    {free(buffer -> pidx); buffer-> pidx  = NULL;}
    if (buffer -> readBuf) {free(buffer -> readBuf);  buffer -> readBuf = NULL;}
    if (buffer -> grids)   {for (int i = 0; i < buffer->terms; i++){if (buffer -> grids[i]){fftwf_free(buffer->grids[i]); buffer->grids[i] = NULL;}}}
    if (buffer -> grids)   {free(buffer -> grids);  buffer -> grids = NULL;}
    if (buffer -> peaks)   {free(buffer -> peaks);  buffer -> peaks = NULL;}

    sdsfree(buffer->outBuf); sdsfree(buffer->spectrum);

    for (int i = 0; i < buffer->terms; i++){if (buffer -> gidx && buffer -> gidx[i]) {free(buffer->gidx[i]); buffer->gidx[i] = NULL;}}
    if (buffer->gidx){free(buffer->gidx); buffer-> gidx = NULL;}

    for (int i = 0; i < buffer->terms; i++){if (buffer -> gdist && buffer -> gdist[i]) {free(buffer->gdist[i]); buffer->gdist[i] = NULL;}}
    if (buffer->gdist){free(buffer->gdist); buffer-> gdist = NULL;}
    #ifndef __AVX__
    for (int i = 0; i < buffer->terms; i++){if (buffer -> weights && buffer -> weights[i]) {free(buffer->weights[i]); buffer->weights[i] = NULL;}}
    if (buffer->weights){free(buffer->weights); buffer-> weights = NULL;}
    #endif
}

static inline int alloc_buffer(buffer_t* buffer, int terms, int n, int size, uint32_t gridLen, int npeaks) {
    buffer->memBlockSize = (gridLen + 16) * sizeof(fftwf_complex);
    buffer->len = n; buffer->allocated = true; buffer->terms = terms;
    buffer->spectrum = sdsempty(); buffer->outBuf = sdsempty();
    if (!buffer->x) {buffer->x = aligned_alloc(64, round_buffer(n * sizeof(double)));}
        if (!buffer->x) goto error;
    if (!buffer->y) {buffer->y = aligned_alloc(64, round_buffer(n * sizeof(float)));}
        if (!buffer->y) goto error;
    if (!buffer->dy) {buffer->dy = aligned_alloc(64, round_buffer(n * sizeof(float)));}
        if (!buffer->dy) goto error;
    if (!buffer->pidx) {buffer->pidx = (uint16_t*) malloc(n * sizeof(uint16_t));}
        if (!buffer->pidx) goto error;
    if (!buffer->readBuf) {buffer->readBuf = aligned_alloc(64, round_buffer(size));}
        if (!buffer->readBuf) goto error;
    if (!buffer->grids) {buffer->grids = calloc(terms, sizeof(complex float *));}
        for (int i = 0; i < buffer->terms; i++){buffer->grids[i] = (fftwf_complex*) fftwf_malloc(buffer->memBlockSize);}
        if (!buffer->grids) goto error;
        else buffer->nGrids = terms;
    if (!buffer->peaks){buffer->peaks = calloc(npeaks, sizeof(peak_t));} buffer->nPeaks = 0;
        if (!buffer->peaks) goto error;
    if (!buffer->gidx) {buffer->gidx = calloc(terms, sizeof(uint32_t **));}//
        if (!buffer->gidx) goto error;
        for (int i = 0; i < buffer->terms; i++){
            buffer->gidx[i] = aligned_alloc(64, round_buffer(n * sizeof(uint32_t)));
            if (!buffer->gidx[i]) goto error;
        }

    if (!buffer->gdist) {buffer->gdist = calloc(terms, sizeof(float *));}//
        if (!buffer->gdist) goto error;
        for (int i = 0; i < buffer->terms; i++){
            buffer->gdist[i] = aligned_alloc(64, round_buffer(n * sizeof(float)));
            if (!buffer->gdist[i]) goto error;
        }
    #ifndef __AVX__
    if (!buffer->weights) {buffer->weights = calloc(terms, sizeof(float *));}//
        if (!buffer->weights) goto error;
        for (int i = 0; i < buffer->terms; i++){
            buffer->weights[i] = aligned_alloc(64, round_buffer(16 * n * sizeof(float)));
            if (!buffer->weights[i]) goto error;
        }
    #endif
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
    double tmp = 0;

    // Center the measurement times - increases precision of future computation
    for (unsigned int i = 1; i < buffer->n; i++) {
        buffer->x[i] -= buffer->x[0];
    }
    buffer->x[0] = 0;

    // Initialize sums for regression
    double sumx = 0, sumxsq = 0, sumy = 0, sumxy = 0;
    double lin = 0, c = 0;

    for (unsigned int i = 0; i < buffer->n; i++) {
        sumx += buffer->x[i];               // sum of x
        sumxsq += buffer->x[i] * buffer->x[i]; // sum of x^2
        sumy += buffer->y[i];               // sum of y
        sumxy += buffer->x[i] * buffer->y[i]; // sum of x*y
    }

    // Calculate the denominator for the slope and intercept
    double denum = (buffer->n * sumxsq) - (sumx * sumx);

    // Compute the slope (lin) and intercept (c)
    lin = ((buffer->n * sumxy) - (sumx * sumy)) / denum;
    c = ((sumy * sumxsq) - (sumx * sumxy)) / denum;
    buffer->magnitude = c;

    // Adjust the y values based on the regression line
    for (unsigned int i = 0; i < buffer->n; i++) {
        buffer->y[i] -= (lin * buffer->x[i]) + c;
    }
}

static inline void linregw_buffer(buffer_t* buffer){
    double tmp = 0;

    // Center the measurement times - increases precision of future computation
    for (unsigned int i = 1; i < buffer->n; i++) {buffer->x[i] -= buffer->x[0];}
    buffer->x[0] = 0;

    // Initialize sums for weighted regression
    double sumx = 0, sumxsq = 0, sumy = 0, sumxy = 0, sumw = 0; // Add sum of weights
    double lin = 0, c = 0, w;

    for (unsigned int i = 0; i < buffer->n; i++) {
        w = 1 / (buffer->dy[i] * buffer->dy[i]); // weight as the inverse of variance
        if (w > 0){
            sumw += w;               // accumulate total weights
            sumx += buffer->x[i] * w;       // weighted sum of x
            sumxsq += buffer->x[i] * buffer->x[i] * w; // weighted sum of x^2
            sumy += buffer->y[i] * w;       // weighted sum of y
            sumxy += (buffer->x[i] * buffer->y[i]) * w; // weighted sum of x*y
        }
        else {buffer->dy[i] = 999.9;}
    }

    // Calculate the denominator for the slope and intercept
    double denum = (sumw * sumxsq) - (sumx * sumx);

    // Compute the slope (lin) and intercept (c)
    lin = ((sumw * sumxy) - (sumx * sumy)) / denum;
    c = ((sumy * sumxsq) - (sumx * sumxy)) / denum;
    buffer->magnitude = c;

    // Adjust the y values based on the regression line
    for (unsigned int i = 0; i < buffer->n; i++) {buffer->y[i] -= (lin * buffer->x[i]) + c;}
}

static inline void preprocess_buffer(buffer_t* buffer, double epsilon){
    linregw_buffer(buffer);

    //adjust the weights
    for(uint32_t i = 0; i < buffer->n; i++){
        buffer->dy[i] *= buffer->dy[i];
        buffer->dy[i] += epsilon;
        buffer->dy[i] = 1.0 / buffer->dy[i];
    }

    double wsum = 0; double wsqsum = 0; double neff = 0;

    for(uint32_t i = 0; i < buffer->n; i++){wsum += fabs(buffer->dy[i]); wsqsum += buffer->dy[i] * buffer->dy[i];}
    neff = ((wsum * wsum) / wsqsum) - 2.0; wsum = 0; // -2 for linear regression
    float nEffInv = 1 / neff;
    //printf("%f\n", neff); // ok
    for(uint32_t i = 0; i < buffer->n; i++){buffer->dy[i] *= buffer->y[i]; wsum += fabs(buffer->dy[i]);} // ok
    wsum = sqrt(neff) / wsum; //sqrt because of square later
    for(uint32_t i = 0; i < buffer->n; i++){buffer->dy[i] *= wsum;} //correct result
}

static inline void read_dat(const char* in_file, buffer_t* buffer) {
    // Open the file
    int fd = open(in_file, O_RDONLY);
    if (fd == -1) {
        perror("Error opening file");
        //return;
    }

    // Get the size of the file
    struct stat file_stat;
    if (fstat(fd, &file_stat) == -1) {
        perror("Failed to get file size");
        close(fd);
        //return;
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
    close(fd); dataBuffer[bytes_read] = '\0';

    // Parse the data
    char* it = dataBuffer;
    char* end = dataBuffer + bytes_read;
    double tempX;
    float tempY, tempDY;
    size_t idx = 0;

    while (it < end && idx < buffer->len) {
        // Parse temporary variables
        tempX = fast_strtod(it, &it); if (it == NULL || it >= end) break; it++;
        tempY = fast_strtof(it, &it); if (it == NULL || it >= end) break; it++;
        tempDY = fast_strtof(it, &it); if (it == NULL || it >= end) break; it++;

        //printf("%f\t%f\t%f\n", tempX, tempY, tempDY); // ok
        buffer->x[idx] = tempX;
        buffer->y[idx] = tempY;
        buffer->dy[idx] = tempDY;
        idx++;
    }

    buffer->n = idx;
}

static inline void append_peak(buffer_t *buff, const int maxPeaks, const double freq, const float magnitude) {
    peak_t appended = {0}; peak_t tmp;
    appended.freq = freq; appended.p = magnitude;
    int idx = buff->nPeaks;

    while (idx > 0 && magnitude > buff->peaks[idx - 1].p) {idx--;}
    if (buff->nPeaks < maxPeaks) {buff->nPeaks++;}
    for (int i = buff->nPeaks - 1; i > idx; i--) {buff->peaks[i] = buff->peaks[i - 1];}
    if (idx < maxPeaks) {buff->peaks[idx] = appended;}
}


#endif
