#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h> // For C99 complex numbers

#include "readout.h"
#include "params.h"

#include "../include/klib/kstring.h"
#include "/usr/local/include/fftw3.h"  // Include FFTW3 header

uint32_t bitCeil(uint32_t n) {
    int exp;
    double mantissa = frexp((double)n, &exp);

    // If the number is a power of two, return it as is
    if (mantissa == 0.5) {return n;}
    // Otherwise, return 2^exp
    return 1 << exp;
}

//wrong results for fmin > 0, grids <= 32768 (2^15) and >= 524288 (2^19). To be fixed (muFFT's bug?)
void process_target(char* in_file, buffer_t* buffer, parameters* params){
    read_dat(kv_A(params->targets, 0).path, buffer); linreg_buffer(buffer); //read the data from .dat file
    int n = 2 + (int)(log(buffer->x[buffer->n-1]) * M_LOG10E); //number of significant digits required for the spectrum

    //adjust the weights
    for(uint32_t i = 0; i < buffer->n; i++){
        buffer->dy[i] *= buffer->dy[i];
        buffer->dy[i] += params->epsilon;
        buffer->dy[i] = 1.0 / buffer->dy[i];
    }

    double wsum = 0; double wsqsum = 0; double neff = 0;

    for(uint32_t i = 0; i < buffer->n; i++){wsum += abs(buffer->dy[i]); wsqsum += buffer->dy[i] * buffer->dy[i];}
    neff = ((wsum * wsum) / wsqsum) - 2.0; wsum = 0; // -2 for linear regression
    //printf("%f\n", neff); // ok
    for(uint32_t i = 0; i < buffer->n; i++){buffer->dy[i] *= buffer->y[i]; wsum += abs(buffer->dy[i]);} // ok
    wsum = sqrt(neff) / wsum; //sqrt because of square later
    for(uint32_t i = 0; i < buffer->n; i++){buffer->dy[i] *= wsum;} //correct result

    double fmax = params -> fmax; double fmin = params -> fmin;
    double fmid = (fmax + fmin) * 0.5; // used to compute the beginning of FFT grid
    double fspan = (fmax - fmin) * (32.0 / 21.0); // used to compute the scale of FFT grid
    uint32_t nsteps = (uint32_t)((double)params->nterms * (double)params->oversamplingFactor * fspan * buffer->x[buffer->n - 1] * 0.5);
    uint32_t gridLen = bitCeil(nsteps);

    printf("Grid size to be allocated: %i\n", gridLen);

    uint32_t memBlockSize = (gridLen + 16) * sizeof(fftwf_complex);

    buffer->grids[0] = (fftwf_complex*) fftwf_malloc(memBlockSize);
    memset(buffer->grids[0], 0, memBlockSize);

    fftwf_plan plan = fftwf_plan_dft_1d(gridLen, buffer->grids[0], buffer->grids[0], FFTW_FORWARD, FFTW_ESTIMATE);

    for(uint32_t i = 0; i < buffer->n; i++){
        double idx_double = buffer->x[i] * fspan; // ok
        uint32_t idx = (uint32_t)(idx_double); // ok
        double idx_frac = idx_double - (double)idx; // ok
        idx = idx % gridLen; // ok
        fftwf_complex val = (cos(fmid * 2.0 * M_PI * buffer->x[i]) * buffer->dy[i]) + (-sin(fmid * 2.0 * M_PI * buffer->x[i]) * buffer->dy[i]) * I; // ok
        if (idx_frac > 0.01){
            float dst = -7.0 - idx_frac; // ok
            for(uint32_t j = 0; j < 16; j++){
                buffer->grids[0][idx + j] += val * sin(dst * M_PI) * sin(dst * M_PI / 3) / (dst * dst * M_PI * M_PI / 3); //sinc(x) * sinc(x/3)
                dst += 1.0;
            }
        } else {
            buffer->grids[0][idx] += val;
        }
    }
    for(uint32_t j = 0; j < 16; j++){
        buffer->grids[0][j] += buffer->grids[0][j+gridLen];
    }

    fftwf_execute(plan);

    if (params->spectrum){

        //negative half
        for(uint32_t i = gridLen * 43 / 64; i < gridLen; i++)
        {ksprintf(&buffer->spectrum, "%.*f\t%.2f\n", n, fmin + (((double)fmax * (double)((i - (gridLen * 43 / 64)) * 32)) / (double)(gridLen * 21)),
            (crealf(buffer->grids[0][i]) * crealf(buffer->grids[0][i]) + cimagf(buffer->grids[0][i]) * cimagf(buffer->grids[0][i])));}

        //positive half
        for(uint32_t i = 0; i <= gridLen * 21 / 64; i++) // higher half
        {ksprintf(&buffer->spectrum, "%.*f\t%.2f\n", n, fmid + (((double)fmax * (double)((i) * 32)) / (double)(gridLen * 21)),
            (crealf(buffer->grids[0][i]) * crealf(buffer->grids[0][i]) + cimagf(buffer->grids[0][i]) * cimagf(buffer->grids[0][i])));}
    }

    if (params->spectrum){
        // Prepare the output file
        char out_file[256];
        strncpy(out_file, in_file, sizeof(out_file) - 1);
        out_file[sizeof(out_file) - 1] = '\0'; // Ensure null-termination
        char* dot = strrchr(out_file, '.'); // Find the last dot in the file name
        if (dot) {*dot = '\0';}  // Remove the existing extension
        strcat(out_file, ".tsv"); // Append the new extension

        // Write the formatted string to the new .tsv file
        FILE* fp = fopen(out_file, "w");
        if (fp == NULL) {perror("Failed to open file for writing"); return;}

        fprintf(fp, "%s\n", buffer->spectrum.s); fclose(fp); // Write the spectrum string to the file
    }

    if(buffer->grids[0]){fftwf_free(buffer->grids[0]); buffer->grids[0] = NULL;}
    fftwf_destroy_plan(plan);

    if(buffer->spectrum.s){free(buffer->spectrum.s);}
}

#endif
