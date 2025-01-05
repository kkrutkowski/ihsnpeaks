#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readout.h"
#include "params.h"

#include "../include/klib/kstring.h"
#include "../include/mufft/mufft.x86.h"

uint32_t bitCeil(uint32_t n) {
    int exp;
    double mantissa = frexp((double)n, &exp);

    // If the number is a power of two, return it as is
    if (mantissa == 0.5) {return n;}
    // Otherwise, return 2^exp
    return 1 << exp;
}

void process_target(char* in_file, buffer_t* buffer, parameters* params){
    read_dat(kv_A(params->targets, 0).path, buffer); linreg_buffer(buffer); //read the data from .dat file
    int n = 2 + (int)(log(buffer->x[buffer->n-1]) * M_LOG10E); //number of significant digits required for the spectrum

        // double number = 0.123456789; printf("\n%i\n", n); // placeholder line
        // ksprintf(&buffer->spectrum, "%.*f", n, number); printf("Formatted string: %s\n", buffer->spectrum.s); //prinft test

    double fmax = params -> fmax; double fmin = params -> fmin;
    double fmid = (fmax + fmin) * 0.5; // used to compute the beginning of FFT grid
    double fspan = (fmax - fmin) * 32 / 21; // used to compute the scale of FFT grid
    uint32_t nsteps = (uint32_t)((double)params->nterms * (double)params->oversamplingFactor * fspan * buffer->x[buffer->n - 1] * 0.5);
    nsteps = bitCeil(nsteps);

        printf("Number of target frequencies: %i\n", nsteps);
    uint32_t gridSize = nsteps;
        printf("Grid size to be allocated: %i\n", gridSize);

    buffer->grids[0] = mufft_alloc((nsteps + 16) * sizeof(complex float));
    printf("%i\n", (nsteps + 16) * sizeof(complex float));
    buffer->plan = mufft_create_plan_1d_c2c(nsteps, MUFFT_FORWARD, 0);
    mufft_execute_plan_1d(buffer->plan, buffer->grids[0], buffer->grids[0]);


    if(buffer->grids[0]){free(buffer->grids[0]);}
    mufft_free_plan_1d(buffer->plan);



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

    if(buffer->spectrum.s){free(buffer->spectrum.s);}
}



#endif
