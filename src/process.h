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

#include "../include/fast_convert.h"
#include "../include/klib/kstring.h"
#include "/usr/local/include/fftw3.h"  // Include FFTW3 header

unsigned int custom_ftoa(float v, int size, char *line) {
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1) {
        line[0] = '0'; line[1] = '.';
        for (int i = 0; i < size; i++) {line[2 + i] = '0';}
        line[4] = '\0'; return 4; // Return the length of the string
    } else {return fast_ftoa(v, new_size, line);}
}

unsigned int custom_dtoa(double v, int size, char *line){
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1){return fast_dtoa(0.0, new_size, line);}
    else {return fast_dtoa(v, new_size, line);}
}

//wrong results for fmin > 0, grids <= 32768 (2^15) and >= 524288 (2^19). To be fixed (muFFT's bug?)
void process_target(char* in_file, buffer_t* buffer, parameters* params){
    read_dat(kv_A(params->targets, 0).path, buffer); linreg_buffer(buffer); //read the data from .dat file
    int n = 1 + (int)(log10(buffer->x[buffer->n-1] * (double)(params->oversamplingFactor * params->nterms))); //number of significant digits required for the spectrum

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

    double fmin = params->fmin; double fstep = 1.0 / (double)(params->nterms * (double)params->oversamplingFactor * buffer->x[buffer->n - 1] * 0.5);
    double fspan = (double)(params->gridLen) * fstep; // used to compute the scale of FFT grid
    double fjump = fspan * (21.0/32.0); //used to switch to next transform
    double fmax = fmin + fjump;
    double fmid = (fmax + fmin) * 0.5; // used to compute the beginning of FFT grid
    uint32_t nsteps = (uint32_t)(params->gridLen);
    uint32_t gridLen = params->gridLen;

    printf("Grid length: %i\n", params->gridLen);
    printf("Number of target frequencies: %i\n", (int)((params->fmax - params->fmin)/fstep));

    uint32_t memBlockSize = (params->gridLen + 16) * sizeof(fftwf_complex);

    buffer->grids[0] = (fftwf_complex*) fftwf_malloc(memBlockSize);

    //fftwf_plan plan = fftwf_plan_dft_1d(gridLen, buffer->grids[0], buffer->grids[0], FFTW_FORWARD, FFTW_MEASURE);
    fftwf_plan plan = fftwf_plan_dft_1d(gridLen, buffer->grids[0], buffer->grids[0], FFTW_FORWARD, FFTW_ESTIMATE);

    // Single buffer to hold the converted strings
    char stringBuff[32];  // Adjust size as needed
    double freq = 0; float magnitude = 0;

    double invGridLen = 1.0 / (double)(gridLen);
    uint32_t shift = (gridLen * 43 / 64);

    while(fmin < params->fmax){

        memset(buffer->grids[0], 0, memBlockSize);

        for(uint32_t i = 0; i < buffer->n; i++){
            double idx_double = buffer->x[i] * fspan; // ok
            uint32_t idx = (uint32_t)(idx_double); // ok
            double idx_frac = idx_double - (double)idx; // ok
            idx = idx % gridLen; // ok
            fftwf_complex val = cexp(-2.0 * I * M_PI * fmid * buffer->x[i]) * buffer->dy[i]; // ok
            if (idx_frac > 0.01){
                float dst = -7.0 - idx_frac; // ok
                for(uint32_t j = 0; j < 16; j++){//replace with vectorized sinc
                    buffer->grids[0][idx + j] += val * sin(dst * M_PI) * sin(dst * M_PI / 3) / (dst * dst * M_PI * M_PI / 3); //sinc(x) * sinc(x/3)
                    dst += 1.0;
                }
            } else {
                buffer->grids[0][idx] += val;
            }
        }
    for(uint32_t j = 0; j < 16; j++){buffer->grids[0][j] += buffer->grids[0][j+gridLen];}

    fftwf_execute(plan);

        for (uint32_t i = shift + 1; i < gridLen; i++) { // negative half
            // Convert the first column value using dtoa
            freq = fmin + ((double)((i) - shift) * invGridLen * fspan); if (freq > params->fmax){goto end;}
            magnitude = crealf(buffer->grids[0][i]) * crealf(buffer->grids[0][i]) + cimagf(buffer->grids[0][i]) * cimagf(buffer->grids[0][i]);
            if (params->spectrum){
                custom_dtoa(freq, n, stringBuff);
                // Append the formatted string to the kstring buffer
                kputs(stringBuff, &buffer->spectrum);
                kputc('\t', &buffer->spectrum);
                // Convert the second column value using ftoa
                custom_ftoa(magnitude, 2, stringBuff);  // 2 decimal places
                // Append the formatted string to the kstring buffer
                kputs(stringBuff, &buffer->spectrum);
                kputc('\n', &buffer->spectrum);
            }
        }

        for (uint32_t i = 0; i <= gridLen * 21 / 64; i++) { // positive half
            // Convert the first column value using dtoa
            freq = fmid + ((double)(i) * invGridLen * fspan);  if (freq > params->fmax){goto end;}
            magnitude = crealf(buffer->grids[0][i]) * crealf(buffer->grids[0][i]) + cimagf(buffer->grids[0][i]) * cimagf(buffer->grids[0][i]);
            if (params->spectrum){
                custom_dtoa(freq, n, stringBuff);
                // Append the formatted string to the kstring buffer
                kputs(stringBuff, &buffer->spectrum);
                kputc('\t', &buffer->spectrum);
                // Convert the second column value using ftoa
                custom_ftoa(magnitude, 2 , stringBuff);  // 2 decimal places
                // Append the formatted string to the kstring buffer
                kputs(stringBuff, &buffer->spectrum);
                kputc('\n', &buffer->spectrum);
            }
        }
        fmin += fjump; fmid += fjump; fmax += fjump;
    } end:


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
