#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <complex.h> // For C99 complex numbers

#include "utils/simd.h"
#include "utils/readout.h"
#include "params.h"

#include "../include/sds.h"
#include "../include/fast_convert.h"
#include "/usr/local/include/fftw3.h"  // Include FFTW3 header

inline float correctPower(float K, float nInv) {
    float term1 = ((2.0 * K) - (K * K)) * (0.25 * nInv);
    float term2 = ((24.0 * K) - (132.0 * K * K) + (76.0 * K * K * K) - (9.0 * K * K * K * K)) * (nInv * nInv * 3.4722222e-3); // 1/288
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

static inline void appendFreq(double freq, float magnitude, int precision, sds* buffer, char* stringBuff) {
    // Convert frequency to string and append to the buffer
    custom_dtoa(freq, precision, stringBuff);
    *buffer = sdscat(*buffer, stringBuff);
    *buffer = sdscatlen(*buffer, "\t", 1);  // Add the separator

    // Convert magnitude to string and append to the buffer
    custom_ftoa(magnitude, 2, stringBuff);
    *buffer = sdscat(*buffer, stringBuff);
    *buffer = sdscatlen(*buffer, "\n", 1);  // Add a newline character
}

static inline void write_tsv(buffer_t *buffer, char* in_file){
        char out_file[256];
        strncpy(out_file, in_file, sizeof(out_file) - 1);
        out_file[sizeof(out_file) - 1] = '\0'; // Ensure null-termination
        char* dot = strrchr(out_file, '.'); // Find the last dot in the file name
        if (dot) {*dot = '\0';}  // Remove the existing extension
        strcat(out_file, ".tsv"); // Append the new extension

        // Write the formatted string to the new .tsv file
        FILE* fp = fopen(out_file, "w");
        if (fp == NULL) {perror("Failed to open file for writing"); return;}

        fprintf(fp, "%s\n", buffer->spectrum); fclose(fp); // Write the spectrum string to the file
        sdsclear(buffer->spectrum);
};

void print_peaks(buffer_t *buffer, parameters *params, int n, char *stringBuff, char *in_file) {
    int i = 0;
    // Append file information to the output buffer
    buffer->outBuf = sdscat(buffer->outBuf, "File: ");
    buffer->outBuf = sdscat(buffer->outBuf, in_file);
    buffer->outBuf = sdscat(buffer->outBuf, "    NofData: ");
    custom_ftoa(buffer->n, 0, stringBuff);  // Convert integer to string using custom_ftoa
    buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
    buffer->outBuf = sdscat(buffer->outBuf, "    Average: ");
    custom_ftoa(buffer->magnitude, 2, stringBuff);  // Convert float to string using custom_ftoa
    buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
    buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);  // Append a newline

    // Append the header for the peaks table
    buffer->outBuf = sdscat(buffer->outBuf, "   f[1/d]\tlog(p)\tAmp\tChi^2\n");

    // Append each peak's information to the output buffer
    while (i < params->npeaks && buffer->peaks[i].p > 0) {
        // Convert frequency to string using custom_ftoa
        custom_dtoa(buffer->peaks[i].freq, n, stringBuff);
        buffer->outBuf = sdscat(buffer->outBuf, "   ");
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);

        // Convert log(p) to string using custom_ftoa
        custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);

        // Convert amplitude to string using custom_ftoa
        custom_ftoa(buffer->peaks[i].amp, 2, stringBuff);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);

        // Convert chi^2 to string using custom_ftoa
        custom_ftoa(buffer->peaks[i].chi2, 2, stringBuff);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);

        i++;
    }
    printf("%s", buffer->outBuf);
    sdsclear(buffer->outBuf);
}

void fprint_buffer(buffer_t *buffer, parameters *params) {
    pthread_mutex_lock(&params->mutex);
    FILE *file = fopen(params->outFile, "a");  // Open the file in append mode
    if (file == NULL) {pthread_mutex_unlock(&params->mutex); perror("Failed to open file for appending"); return;}
    fprintf(file, "%s", buffer->outBuf); fflush(stdout); fclose(file);
    sdsclear(buffer->outBuf);
    pthread_mutex_unlock(&params->mutex);
}

void append_peaks(buffer_t *buffer, parameters *params, int n, char *stringBuff, char *in_file) {
    int i = 0;
        if (buffer->nPeaks > 0){
        // Append file information to the output buffer
        buffer->outBuf = sdscat(buffer->outBuf, in_file);
        custom_dtoa(buffer->peaks[0].freq, n, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        custom_ftoa(buffer->peaks[0].p * M_LOG10E, 2, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        //buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);

        // Append each peak's information to the output buffer
        while (i < params->npeaks && buffer->peaks[i].p > 0) {
            // Convert frequency to string using custom_ftoa
            buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
            custom_dtoa(buffer->peaks[i].freq, n, stringBuff);
            buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
            buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
            // Convert log(p) to string using custom_ftoa
            custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
            buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
            buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
            i++;
        }
        buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);
        buffer->loc_iter += 1;
        if(buffer->loc_iter == 0){ //append to the file only if the local counter oveflows
            // Append the buffer content to the file
            fprint_buffer(buffer, params); fflush(stdout);
        }
    }
}

//wrong results for fmin > 0, grids <= 32768 (2^15) and >= 524288 (2^19). To be fixed (muFFT's bug?)
void process_target(char* in_file, buffer_t* buffer, parameters* params, const bool batch){
    read_dat(in_file, buffer); linreg_buffer(buffer); //read the data from .dat file
    const int n = 1 + (int)(log10(buffer->x[buffer->n-1] * (double)(params->oversamplingFactor * params->nterms))); //number of significant digits required for the spectrum
    const float treshold = params->treshold * M_LN10; memset(buffer->peaks, 0, params->npeaks * sizeof(peak_t));

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

    if(!batch && params->debug){printf("\tNumber of target frequencies: %i\n\n", (int)((params->fmax - params->fmin)/fstep));}

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
                fftwf_complex val = cexp(-2.0 * M_PI * freqFactor * I * fmid * buffer->x[i]) * buffer->dy[i]; // cexp takes ~ 15% of the compute time
                if (buffer->gdist[t][i] < 0.01){buffer->grids[t][buffer->gidx[t][i]] += val;}
                else if (buffer->gdist[t][i] > 0.99){buffer->grids[t][buffer->gidx[t][i+1]] += val;}
                else {
                #ifdef __AVX__
                VEC weights[2];
                generateWeights(buffer->gdist[t][i], &weights[0], &weights[1]);
                for(uint32_t j = 0; j < 8; j++){
                    buffer->grids[t][buffer->gidx[t][i] + j] += val * weights[0].data[j];
                    buffer->grids[t][buffer->gidx[t][i] + 8 + j] += val * weights[1].data[j];
                }
                #else
                    float dst = -7.0 - buffer->gdist[t][i];
                    for(uint32_t j = 0; j < 16; j++){
                        buffer->grids[t][buffer->gidx[t][i] + j] += val * 3.0f * sinf(dst * M_PI) * sinf(dst * M_PI * 0.33333333f) / (dst * dst * M_PI * M_PI);
                        dst += 1.0;
                    }
                #endif
                }
            }
        for(uint32_t j = 0; j < 16; j++){buffer->grids[t][j] += buffer->grids[t][j+gridLen];}

        fftwf_execute_dft(params->plan, buffer->grids[t], buffer->grids[t]);
    }
        float magnitudes[3] = {0};
        for(int t = 0; t < buffer->terms; t++){magnitudes[0] += sabs(buffer->grids[t][shift]); magnitudes[1] += sabs(buffer->grids[t][shift+1]);}

        for (uint32_t i = shift + 1; i < gridLen; i++) { // negative half
            freq = fmin + ((double)((i) - shift) * invGridLen * fspan); if (freq > params->fmax){goto end;}
            magnitudes[2] = 0;
            for(int t = 0; t < buffer->terms; t++){magnitudes[2] += sabs(buffer->grids[t][i+1]);}
            //magnitudes[1] += correctPower(sabs(buffer->grids[t][i]), nEffInv);
            if (params->spectrum){appendFreq(freq, magnitudes[1], n, &buffer->spectrum, &stringBuff[0]);}
            if (magnitudes[1] > treshold && magnitudes[1] > magnitudes[0] && magnitudes[1] > magnitudes[2]){append_peak(buffer, params->npeaks, freq, magnitudes[1]);}
            magnitudes[0] = magnitudes[1]; magnitudes[1] = magnitudes[2]; //reuse the results in the next iteration
        }

        for(int t = 0; t < buffer->terms; t++){magnitudes[1] += sabs(buffer->grids[t][0]);}
        if (params->debug && params->spectrum) {buffer->spectrum = sdscat(buffer->spectrum, "break\n");}

        for (uint32_t i = 0; i <= gridLen * 21 / 64; i++) { // positive half
            freq = fmid + ((double)(i) * invGridLen * fspan);  if (freq > params->fmax){goto end;}
            magnitudes[2] = 0;
            for(int t = 0; t < buffer->terms; t++){magnitudes[2] += sabs(buffer->grids[t][i+1]);}
            //magnitudes[1] += correctPower(sabs(buffer->grids[t][i]), nEffInv);
            if (params->spectrum){appendFreq(freq, magnitudes[1], n, &buffer->spectrum, &stringBuff[0]);}
            if (magnitudes[1] > treshold && magnitudes[1] > magnitudes[0] && magnitudes[1] > magnitudes[2]){append_peak(buffer, params->npeaks, freq, magnitudes[1]);}
            magnitudes[0] = magnitudes[1]; magnitudes[1] = magnitudes[2]; //reuse the results in the next iteration
        }
        fmin += fjump; fmid += fjump; fmax += fjump; if (params->debug && params->spectrum) {buffer->spectrum = sdscat(buffer->spectrum, "break\n");}
    } end:

    if (!batch) {print_peaks(buffer, params, n, stringBuff, in_file);}
    else {append_peaks(buffer, params, n, stringBuff, in_file);}
    if (params->spectrum){write_tsv(buffer, in_file);}
    buffer->nPeaks = 0; //reset the data for buffer reuse
}



void process_targets(void *data, long i, int thread_id) {
    parameters *params = (parameters *)data;

    int permile; permile = kv_size(params->targets) / 1000; if(permile == 1){permile += 1;}
    pthread_mutex_lock(&params->counter_mutex); params->iter_count += 1; pthread_mutex_unlock(&params->counter_mutex);
    if (params->iter_count % permile == 0 || params->iter_count == kv_size(params->targets) - 1) {
        float progress = (float)(params->iter_count + 1) * 100.0 / (float)(kv_size(params->targets) - 1);
        printf("Computation in progress: %.1f%% complete\r", progress);
        fflush(stdout);
    }
    if (!params->buffers[thread_id]->allocated) {alloc_buffer(params->buffers[thread_id], params->nterms, params->maxLen, params->maxSize, params->gridLen, params->npeaks);
        //printf("Allocating buffer for thread %i\n", thread_id); //ok
    }
    //printf("Processing iteration %ld on thread %d\t%s\n", i, thread_id, kv_A(params->targets, i).path); //ok
    // Process the data for this iteration
    process_target(kv_A(params->targets, i).path, params->buffers[thread_id], params, true);
}


#endif
