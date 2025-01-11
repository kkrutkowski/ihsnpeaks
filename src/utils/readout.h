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

#include "../../include/klib/kstring.h"
#include "../../include/fast_convert.h"

//#include "../include/mufft/mufft.x86.h"
#include "/usr/local/include/fftw3.h"

static inline size_t round_buffer(size_t size) {return (size + 63) & ~63;}

typedef struct {
    double freq;
    float  p;
    float  amp;
    float  chi2;
} peak_t;

typedef struct {
    bool allocated;
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

    kstring_t spectrum;
    kstring_t outBuf;
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

    for (int i = 0; i < buffer->terms; i++){if (buffer -> gidx && buffer -> gidx[i]) {free(buffer->gidx[i]); buffer->gidx[i] = NULL;}}
    if (buffer->gidx){free(buffer->gidx); buffer-> gidx = NULL;}

    for (int i = 0; i < buffer->terms; i++){if (buffer -> gdist && buffer -> gdist[i]) {free(buffer->gdist[i]); buffer->gdist[i] = NULL;}}
    if (buffer->gdist){free(buffer->gdist); buffer-> gdist = NULL;}
}

static inline int alloc_buffer(buffer_t* buffer, int terms, int n, int size, uint32_t gridLen, int npeaks) {
    buffer->memBlockSize = (gridLen + 16) * sizeof(fftwf_complex);
    buffer->len = n; buffer->allocated = true; buffer->terms = terms;
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

    if (!buffer->gdist) {buffer->gdist = calloc(terms, sizeof(float **));}//
        if (!buffer->gdist) goto error;
        for (int i = 0; i < buffer->terms; i++){
            buffer->gdist[i] = aligned_alloc(64, round_buffer(n * sizeof(float)));
            if (!buffer->gdist[i]) goto error;
        }
    return 0;

error:
    free_buffer(buffer);
    fprintf(stderr, "Failed to allocate buffer\n");
    return -1;
}

void print_buffer(buffer_t* buffer) {
    for (int i = 0; i < buffer->n; i++) {
        printf("%.2f\t%.2f\t%.2f\n", buffer->x[i], buffer->y[i], buffer->dy[i]);
    }
}

void linreg_buffer(buffer_t* buffer){
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

void read_dat(const char* in_file, buffer_t* buffer) {
    // Open the file
    int fd = open(in_file, O_RDONLY);
    if (fd == -1) {
        perror("Error opening file");
        return;
    }

    // Get the size of the file
    struct stat file_stat;
    if (fstat(fd, &file_stat) == -1) {
        perror("Failed to get file size");
        close(fd);
        return;
    }
    size_t file_size = file_stat.st_size;

    // Use existing buffer to avoid reallocation
    char* dataBuffer = buffer->readBuf;

    // Read the file
    ssize_t bytes_read = read(fd, dataBuffer, file_size);
    if (bytes_read == -1) {perror("Failed to read file"); close(fd); return;
    } else if (bytes_read == 0) {
        fprintf(stderr, "No data read from file\n"); close(fd); return;
    }

    dataBuffer[bytes_read] = '\0'; // Null-terminate the buffer

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

    // Close the file
    close(fd); buffer->n = idx;
}


#endif
