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

#include "../include/klib/kvec.h"
#include "../include/fast_convert.h"
#include "../include/mufft/mufft.x86.h"

//struct parameters;

typedef struct {
    bool allocated;
    int        len;
    int          n;

    char*  readBuf;
    double*      x;
    float*       y;
    float*      dy;
    uint32_t* gidx; //grid indices for NFFT
    uint16_t* pidx; //phase indices for counting sort

    //D_VEC*    bufd1;
    //D_VEC*    bufd2;

    float magnitude;
    int gridSize;
    complex float ** grids; //used to compute the FFT

    mufft_plan_1d *muplan; //FFT plan // = mufft_create_plan_1d_c2c(N, MUFFT_FORWARD, flags);
    // https://github.com/Themaister/muFFT/blob/master/bench.c
} buffer_t;

static inline void free_buffer (buffer_t* buffer) {
    if (buffer -> x)     {free(buffer -> x);  buffer -> x = NULL;}
    if (buffer -> y)     {free(buffer -> y);  buffer -> y = NULL;}
    if (buffer -> dy)    {free(buffer -> dy); buffer-> dy  = NULL;}
    if (buffer -> gidx)  {free(buffer -> gidx); buffer-> gidx  = NULL;}
    if (buffer -> pidx)  {free(buffer -> pidx); buffer-> pidx  = NULL;}
    if (buffer -> readBuf) {free(buffer -> readBuf);  buffer -> readBuf = NULL;}
}

static inline size_t round_buffer(size_t size) {return (size + 63) & ~63;}

static inline int alloc_buffer(buffer_t* buffer, int n, int size) {
    buffer->len = n; buffer->allocated = true;
    if (!buffer->x) {buffer->x = aligned_alloc(64, round_buffer(n * sizeof(double)));}
        if (!buffer->x) goto error;
    if (!buffer->y) {buffer->y = aligned_alloc(64, round_buffer(n * sizeof(float)));}
        if (!buffer->y) goto error;
    if (!buffer->dy) {buffer->dy = aligned_alloc(64, round_buffer(n * sizeof(float)));}
        if (!buffer->dy) goto error;
    if (!buffer->gidx) {buffer->gidx = aligned_alloc(64, round_buffer(n * sizeof(uint32_t)));}
        if (!buffer->gidx) goto error;
    if (!buffer->pidx) {buffer->pidx = (uint16_t*) malloc(n * sizeof(uint16_t));}
        if (!buffer->pidx) goto error;
    if (!buffer->readBuf) {buffer->readBuf = aligned_alloc(64, round_buffer(size));}
        if (!buffer->readBuf) goto error;
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

    buffer->magnitude = lin;

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
