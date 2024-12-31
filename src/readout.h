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

#include "../include/klib/kvec.h"
#include "../include/fast_convert.h"

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

    //S_VEC*    bufs1;
    //S_VEC*    bufs2;
    //D_VEC*    bufd1;
    //D_VEC*    bufd2;

    float ** grids; //used to compute the FFT
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
        // Parse tempX
        tempX = fast_strtod(it, &it);
        if (it == NULL || it >= end) break; it++;
        tempY = fast_strtof(it, &it); if (it == NULL || it >= end) break; it++;
        tempDY = fast_strtof(it, &it); if (it == NULL || it >= end) break; it++;

        buffer->x[idx] = tempX;
        buffer->y[idx] = tempY;
        buffer->dy[idx] = tempDY;
        idx++; buffer->n = idx;
    }

    // Close the file
    close(fd);
}




#endif
