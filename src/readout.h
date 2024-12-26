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

typedef struct {
    double*      x;
    float*       y;
    float*      dy;
    uint32_t* gidx; //grid indices for NFFT
    uint16_t* pidx; //phase indices for counting sort
    //S_VEC*    bufs1;
    //S_VEC*    bufs2;
    //D_VEC*    bufd1;
    //D_VEC*    bufd2;
} buffer_t;

inline int alloc_buffer (buffer_t* buffer, int n) {
    if (!buffer -> x)     {buffer -> x = aligned_alloc(64, n * sizeof(double));
        if (!buffer->x) return fprintf(stderr, "Failed to allocate buffer\n"), -1;}
    if (!buffer -> y)     {buffer -> y = aligned_alloc(64, n * sizeof(float));
        if (!buffer->y) return fprintf(stderr, "Failed to allocate buffer\n"), -1;}
    if (!buffer -> dy)    {buffer -> dy = aligned_alloc(64, n * sizeof(float));
        if (!buffer->dy) return fprintf(stderr, "Failed to allocate buffer\n"), -1;}
    if (!buffer -> gidx)  {buffer -> gidx = aligned_alloc(64, n * sizeof(uint32_t));
        if (!buffer->gidx) return fprintf(stderr, "Failed to allocate buffer\n"), -1;}
    if (!buffer -> pidx)  {buffer -> pidx = aligned_alloc(64, 1024 * sizeof(uint32_t));
        if (!buffer->pidx) return fprintf(stderr, "Failed to allocate buffer\n"), -1;}
return 0;}

inline void free_buffer (buffer_t* buffer) {
    if (buffer -> x)     {free(buffer -> x);  buffer -> x = NULL;}
    if (buffer -> y)     {free(buffer -> y);  buffer -> y = NULL;}
    if (buffer -> dy)    {free(buffer -> dy); buffer-> dy  = NULL;}
    if (buffer -> gidx)  {free(buffer -> gidx); buffer-> gidx  = NULL;}
    if (buffer -> pidx)  {free(buffer -> pidx); buffer-> pidx  = NULL;}
}

typedef struct {
    kvec_t(double) x;
    kvec_t(float)  y;
    kvec_t(float) dy;
} star_t;



inline void read_dat(const char* in_file, buffer_t* buffer) {
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

    // Allocate a buffer to read the file
    char* dataBuffer = (char*)malloc(file_size + 1); // +1 for null-termination
    if (!dataBuffer) {
        fprintf(stderr, "Memory allocation failed\n");
        close(fd);
        return;
    }

    // Read the file
    ssize_t bytes_read = read(fd, dataBuffer, file_size);
    if (bytes_read == -1) {
        perror("Failed to read file");
        free(dataBuffer);
        close(fd);
        return;
    } else if (bytes_read == 0) {
        fprintf(stderr, "No data read from file\n");
        free(dataBuffer);
        close(fd);
        return;
    }

    dataBuffer[bytes_read] = '\0'; // Null-terminate the buffer

    char* it = dataBuffer;
    char* end = dataBuffer + bytes_read;
    double tempX;
    float tempY, tempDY;
    size_t idx = 0;

    while (it < end) {
        // Parse tempX
        tempX = fast_strtod(it, &it);
        if (it == NULL || it >= end) break; it++;
        tempY = fast_strtof(it, &it); if (it == NULL || it >= end) break; it++;
        tempDY = fast_strtof(it, &it); if (it == NULL || it >= end) break; it++;

        if (!buffer->x || !buffer->y || !buffer->dy) {
            fprintf(stderr, "Memory reallocation failed\n");
            free(dataBuffer);
            close(fd);
            return;
        }

        buffer->x[idx] = tempX;
        buffer->y[idx] = tempY;
        buffer->dy[idx] = tempDY;
        idx++;
    }

    // Free temporary buffer
    free(dataBuffer);

    // Close the file
    close(fd);
}




#endif
