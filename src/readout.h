#ifndef READOUT_H
#define READOUT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
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

typedef struct {
    kvec_t(double) x;
    kvec_t(float)  y;
    kvec_t(float) dy;
} star_t;










#endif
