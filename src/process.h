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

float correctPower(float K, float nInv) {
    float term1 = ((2.0 * K) - (K * K)) * (0.25 * nInv);
    float term2 = ((24.0 * K) - (132.0 * K * K) + (76.0 * K * K * K) - (9.0 * K * K * K * K)) * (nInv * nInv / 288.0);
    float inside_log = logf(1 + term1 - term2);
return K - inside_log;}

inline float sabs(complex float z) {return (creal(z) * creal(z)) + (cimag(z) * cimag(z));}

inline unsigned int custom_ftoa(float v, int size, char *line) {
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1) {
        line[0] = '0'; line[1] = '.';
        for (int i = 0; i < size; i++) {line[2 + i] = '0';}
        line[4] = '\0'; return 4; // Return the length of the string
    } else {return fast_ftoa(v, new_size, line);}
}

inline unsigned int custom_dtoa(double v, int size, char *line){
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1){return fast_dtoa(0.0, new_size, line);}
    else {return fast_dtoa(v, new_size, line);}
}

static inline void appendFreq(double freq, float magnitude, int precision, kstring_t* buffer, char* stringBuff) {
    // Convert frequency to string and append to the buffer
    custom_dtoa(freq, precision, stringBuff);
    kputs(stringBuff, buffer);
    kputc('\t', buffer); // Add the separator
    // Convert magnitude to string and append to the buffer
    custom_ftoa(magnitude, 2, stringBuff);
    kputs(stringBuff, buffer);
    kputc('\n', buffer); // Add a newline character
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
    float nEffInv = 1 / neff;
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

    // Single buffer to hold the converted strings
    char stringBuff[32];  // Adjust size as needed

    double invGridLen = 1.0 / (double)(gridLen);
    uint32_t shift = (gridLen * 43 / 64);

    for(int t = 0; t < buffer->terms; t++){
        for(uint32_t i = 0; i < buffer->n; i++){
            double freqFactor = (double)(t+1);
            double idx_double = buffer->x[i] * fspan * freqFactor;
            uint32_t idx_tmp =  (uint32_t)(idx_double);
            buffer->gdist[t][i] = idx_double - (double)idx_tmp;
            buffer->gidx[t][i] = idx_tmp % gridLen;
        }
    }

    while(fmin < params->fmax){
        double freq = 0;
        for(int t = 0; t < buffer->terms; t++){ //compute the NFFTs required for the data
            memset(buffer->grids[t], 0, buffer->memBlockSize);
            double freqFactor = (double)(t+1); //printf("%.2f\n", freqFactor);
            for(uint32_t i = 0; i < buffer->n; i++){
                fftwf_complex val = cexp(-2.0 * M_PI * freqFactor * I * fmid * buffer->x[i]) * buffer->dy[i]; // ok
                if (buffer->gdist[t][i] > 0.01){
                    float dst = -7.0 - buffer->gdist[t][i]; // ok
                    for(uint32_t j = 0; j < 16; j++){//replace with vectorized sinc
                        buffer->grids[t][buffer->gidx[t][i] + j] += val * 3.0f * sinf(dst * M_PI) * sinf(dst * M_PI * 0.33333333f) / (dst * dst * M_PI * M_PI); //sinc(x) * sinc(x/3)
                        dst += 1.0;
                    }
                } else {
                    buffer->grids[t][buffer->gidx[t][i]] += val;
                }
            }
        for(uint32_t j = 0; j < 16; j++){buffer->grids[t][j] += buffer->grids[t][j+gridLen];}

        fftwf_execute_dft(params->plan, buffer->grids[t], buffer->grids[t]);
    }

        for (uint32_t i = shift + 1; i < gridLen; i++) { // negative half
            freq = fmin + ((double)((i) - shift) * invGridLen * fspan); if (freq > params->fmax){goto end;}
             float magnitude = 0;
            for(int t = 0; t < buffer->terms; t++){magnitude += correctPower(sabs(buffer->grids[t][i]), nEffInv);}
            if (params->spectrum){appendFreq(freq, magnitude, n, &buffer->spectrum, &stringBuff[0]);}
        }

        for (uint32_t i = 0; i <= gridLen * 21 / 64; i++) { // positive half
            freq = fmid + ((double)(i) * invGridLen * fspan);  if (freq > params->fmax){goto end;}
            float magnitude = 0;
            for(int t = 0; t < buffer->terms; t++){magnitude += correctPower(sabs(buffer->grids[t][i]), nEffInv);}
            if (params->spectrum){appendFreq(freq, magnitude, n, &buffer->spectrum, &stringBuff[0]);}
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

    if(buffer->spectrum.s){free(buffer->spectrum.s);}
}

#endif
