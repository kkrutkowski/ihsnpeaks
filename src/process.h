#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <complex.h> // For C99 complex numbers

#include "params.h"

#include "utils/simd.h"
#include "utils/readout.h"
#include "utils/strings.h"
//#include "utils/convolution.h"


#include <fftw3.h>  // Include FFTW3 header


static inline float correctPower(float K, float nInv) {
    float term1 = ((2.0 * K) - (K * K)) * (0.25 * nInv);
    float term2 = ((24.0 * K) - (132.0 * K * K) + (76.0 * K * K * K) - (9.0 * K * K * K * K)) * (nInv * nInv * 3.4722222e-3); // 1/288
    float inside_log = 1 + term1 - term2;
    VEC tmp;
    tmp.data[0] = inside_log;
    inside_log = ln_ps(tmp).data[0];
return K - inside_log;}

inline float sabs(complex float z) {return (creal(z) * creal(z)) + (cimag(z) * cimag(z));}

void process_target(char* in_file, buffer_t* buffer, parameters* params, const bool batch){
    read_dat(in_file, buffer); preprocess_buffer(buffer, params->epsilon, params->mode); //read the data from .dat file

    /* //doesn't work for some reason - using direct indexing instead
    if(params->mode > 0){ //prepare the values array for the F-test
        volatile kvpair* tmp = (volatile kvpair*)(buffer->buf[0]);
        for(int i = 0; i < buffer->n; i++){tmp[i].parts.val = buffer->y[i];}
    }*/

    int prewhitening_iter = 0;

    //a temporary fix of the 'disappearing' r/p values - to be replaced
    double *backup_r = calloc(params->npeaks, sizeof(double)); double *backup_p = calloc(params->npeaks, sizeof(double));

    const int n = 1 + (int)(log10(buffer->x[buffer->n-1] * (double)(params->oversamplingFactor * params->nterms))); //number of significant digits required for the spectrum
    const float threshold = params->threshold;


    double fmin = params->fmin;
    double fstep = 1.0 / (double)(params->nterms * (double)params->oversamplingFactor * buffer->x[buffer->n - 1] * 0.5);
    double fspan = (double)(params->gridLen) * fstep; // used to compute the scale of FFT grid
    double fjump = fspan * (21.0/32.0); //used to switch to next transform
    double fmax = fmin + fjump;
    double fmid = (fmax + fmin) * 0.5; // used to compute the beginning of FFT grid
    //uint32_t nsteps = (uint32_t)(params->gridLen);
    uint32_t gridLen = params->gridLen;

    if(!batch && params->debug){printf("\tNumber of target frequencies: %i\n\n", (int)((params->fmax - params->fmin)/fstep));}

    // Single buffer to hold the converted strings
    char stringBuff[32];  // Adjust size as needed

    double invGridLen = 1.0 / (double)(gridLen);
    double df = 1.0 / buffer->x[buffer->n - 1];
    uint32_t shift = (gridLen * 43 / 64);

    for(int t = 0; t < buffer->terms; t++){
        #ifndef __AVX__
        memset(buffer->weights[t], 0, params->maxLen * 16 * sizeof(float));
        #endif
        for(uint32_t i = 0; i < buffer->n; i++){
            double freqFactor = (double)(t+1);
            double idx_double = buffer->x[i] * fspan * freqFactor;
            uint32_t idx_tmp =  (uint32_t)(idx_double);
            buffer->gdist[t][i] = idx_double - (double)idx_tmp;
            buffer->gidx[t][i] = idx_tmp % gridLen;
            #ifndef __AVX__
            if (buffer->gdist[t][i] < 0.01){buffer->weights[t][16*i] = 1;}
            else if (buffer->gdist[t][i] > 0.99){buffer->weights[t][(16 * i)+1] = 1;}
            else {
                float dst = -7.0 - buffer->gdist[t][i];
                for(uint32_t j = 0; j < 16; j++){
                        buffer->weights[t][(16*i)+j] = 3.0f * sinf(dst * M_PI) * sinf(dst * M_PI * 0.33333333f) / (dst * dst * M_PI * M_PI);
                        dst += 1.0;
                    }
                }
            #endif
        }
    }

    prewhiten:
        if (prewhitening_iter > 0){
            fmin = params->fmin; fmax = fmin + fjump; fmid = (fmax + fmin) * 0.5;
            double prewhitening_freq = buffer->peaks[prewhitening_iter - 1].freq;

            //actual prewhitening happens here
            get_r(buffer, prewhitening_freq, NULL, true);
            //printf("%.3f\n", buffer->peaks[prewhitening_iter - 1].r); //ok

            //temporary fix, to be replaced
            backup_r[prewhitening_iter] = buffer->peaks[prewhitening_iter - 1].r;
            backup_p[prewhitening_iter] = buffer->peaks[prewhitening_iter - 1].p;

            double ysum = 0;
            for (unsigned int i = 0; i < buffer->n; i++){ysum += fabs(buffer->y[i]);}
            float norm = (sqrtf(buffer->neff) / (float)(ysum));
            for (unsigned int i = 0; i < buffer->n; i++){buffer->dy[i] = buffer->y[i] * norm;}

            if (buffer->peaks[prewhitening_iter - 1].r < params->r_threshold) {buffer->nPeaks -= 1; goto next;}
        }

    //memset(&buffer->peaks[prewhitening_iter], 0, params->npeaks * sizeof(peak_t));
    buffer->nPeaks = prewhitening_iter;

    while(fmin < params->fmax){
        double freq = 0;

        //compute the RZTs here
        for(int t = 0; t < buffer->terms; t++){ //compute the NFFTs required for the data
            memset(buffer->grids[t], 0, buffer->memBlockSize);
            double freqFactor = (double)(t+1); //printf("%.2f\n", freqFactor);
            for(uint32_t i = 0; i < buffer->n; i++){
                fftwf_complex val = cexp(-2.0 * M_PI * freqFactor * I * fmid * buffer->x[i]) * buffer->dy[i]; // cexp takes ~ 15% of the compute time
                if (buffer->gdist[t][i] < 0.01){buffer->grids[t][buffer->gidx[t][i]] += val;}
                else if (buffer->gdist[t][i] > 0.99){buffer->grids[t][buffer->gidx[t][i+1]] += val;}
                else {
                #ifdef __AVX512F__
                    VEC weights = generateWeights(buffer->gdist[t][i]);
                    for(uint32_t j = 0; j < 16; j++){buffer->grids[t][buffer->gidx[t][i] + j] += val * weights.data[j];}
                #else
                #ifdef __AVX__
                VEC weights[2];
                generateWeights(buffer->gdist[t][i], &weights[0], &weights[1]);
                for(uint32_t j = 0; j < 8; j++){
                    buffer->grids[t][buffer->gidx[t][i] + j] += val * weights[0].data[j];
                    buffer->grids[t][buffer->gidx[t][i] + 8 + j] += val * weights[1].data[j];
                    }
                #else
                for(uint32_t j = 0; j < 16; j++){buffer->grids[t][buffer->gidx[t][i] + j] += val * buffer->weights[t][(i*16)+j];}
                    #endif // __AVX__
                #endif // __AVX512F__
                }
            }
        for(uint32_t j = 0; j < 16; j++){buffer->grids[t][j] += buffer->grids[t][j+gridLen];}

        if (params->mode < 5){fftwf_execute_dft(params->plan, buffer->grids[t], buffer->grids[t]);}
    }

        float NeffInv = 1.0 / buffer->neff;

        // Inteppret the results - first half
        float magnitudes[3] = {0};
        for(int t = 0; t < buffer->terms; t++) {magnitudes[0] += correctPower(sabs(buffer->grids[t][shift]), NeffInv) ; magnitudes[1] += correctPower(sabs(buffer->grids[t][shift+1]), NeffInv);}


        for (uint32_t i = shift + 1; i < gridLen; i++) { // negative half
            freq = fmin + ((double)((i) - shift) * invGridLen * fspan); if (freq > params->fmax){goto end;}
            magnitudes[2] = 0;
            if(params->mode < 5){for(int t = 0; t < buffer->terms; t++){magnitudes[2] += correctPower(sabs(buffer->grids[t][i+1]), NeffInv);}}
            else {magnitudes[2] = get_z(get_r(buffer, freq + (invGridLen * fspan), NULL, false), buffer->n);}
            if (params->spectrum && params->mode < 5){appendFreq(freq, correct_ihs_res(magnitudes[1], params->nterms), n, &buffer->spectrum, &stringBuff[0]);}
            if (params->spectrum && params->mode > 4){appendFreq(freq, magnitudes[1], n, &buffer->spectrum, &stringBuff[0]);}
            if (magnitudes[1] > threshold && magnitudes[1] > magnitudes[0] && magnitudes[1] > magnitudes[2]){append_peak(buffer, params->npeaks, params->mode, freq, magnitudes[1], df);}
            magnitudes[0] = magnitudes[1]; magnitudes[1] = magnitudes[2]; //reuse the results in the next iteration
        }

        for(int t = 0; t < buffer->terms; t++){magnitudes[1] += correctPower(sabs(buffer->grids[t][0]), NeffInv);}
        if (params->debug && params->spectrum) {buffer->spectrum = sdscat(buffer->spectrum, "break\n");}

        for (uint32_t i = 0; i <= gridLen * 21 / 64; i++) { // positive half
            freq = fmid + ((double)(i) * invGridLen * fspan);  if (freq > params->fmax){goto end;}
            magnitudes[2] = 0;
            if(params->mode < 5){for(int t = 0; t < buffer->terms; t++){magnitudes[2] += correctPower(sabs(buffer->grids[t][i+1]), NeffInv);}}
            else {magnitudes[2] = get_z(get_r(buffer, freq + (invGridLen * fspan), NULL, false), buffer->n);}
            if (params->spectrum && params->mode < 5){appendFreq(freq, correct_ihs_res(magnitudes[1], params->nterms), n, &buffer->spectrum, &stringBuff[0]);}
            if (params->spectrum && params->mode > 4){appendFreq(freq, magnitudes[1], n, &buffer->spectrum, &stringBuff[0]);}
            if (magnitudes[1] > threshold && magnitudes[1] > magnitudes[0] && magnitudes[1] > magnitudes[2]){append_peak(buffer, params->npeaks, params-> mode, freq, magnitudes[1], df);}
            magnitudes[0] = magnitudes[1]; magnitudes[1] = magnitudes[2]; //reuse the results in the next iteration
        }

        //increase frequency
        fmin += fjump; fmid += fjump; fmax += fjump; if (params->debug && params->spectrum) {buffer->spectrum = sdscat(buffer->spectrum, "break\n");}
    } end:


    sortPeaks(&buffer->peaks[prewhitening_iter], buffer->nPeaks, buffer, params->mode, df, params->nterms);

    if (params-> prewhiten && prewhitening_iter < buffer->nPeaks) {
        prewhitening_iter += 1;
        buffer->nPeaks = prewhitening_iter;
            //debugging prints for memory loss - do not remove
            //for (int i = 0; i < buffer->nPeaks; i++){printf("%.3f\n", buffer->peaks[i].r);}
            //printf("\n");
        goto prewhiten
    ;}

    next:
    //a temporary memory loss fix
    if (prewhitening_iter > 1){
        for (int i = 1; i <= prewhitening_iter; i++){
            buffer->peaks[i-1].r = backup_r[i];
            buffer->peaks[i-1].p = backup_p[i];
        };
        buffer->nPeaks = prewhitening_iter;
    }

    if (!batch) {print_peaks(buffer, params, n, stringBuff, in_file, params->mode);}
    else {append_peaks(buffer, params, n, stringBuff, in_file, params->mode);}
    if (params->spectrum){write_tsv(buffer, in_file);}
    free(backup_r); free(backup_p); //temporary "fix", to be removed later
    buffer->nPeaks = 0; //reset the data for buffer reuse
}



void process_targets(void *data, long i, int thread_id) {
    parameters *params = (parameters *)data;

    int permile; permile = kv_size(params->targets) / 1000; if(permile == 0){permile += 1;}
    pthread_mutex_lock(&params->counter_mutex); params->iter_count += 1; pthread_mutex_unlock(&params->counter_mutex);
    if (params->iter_count % permile == 0 || params->iter_count == kv_size(params->targets) - 1) {
        float progress = (float)(params->iter_count + 1) * 100.0 / (float)(kv_size(params->targets) - 1);
        printf("Computation in progress: %.1f%% complete\r", progress);
        fflush(stdout);
    }
    if (!params->buffers[thread_id]->allocated) {alloc_buffer(params->buffers[thread_id], params);
        //printf("Allocating buffer for thread %i\n", thread_id); //ok
    }
    //printf("Processing iteration %ld on thread %d\t%s\n", i, thread_id, kv_A(params->targets, i).path); //ok
    // Process the data for this iteration
    process_target(kv_A(params->targets, i).path, params->buffers[thread_id], params, true);
}


#endif
