#ifndef FTEST_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

#include <fdist.h>

#include "common.h"
#include "simd.h"

// Computes the values of self-correlation of the convolved array elements
static inline double corr(int r){
    double L = 2 * r + 1;
return(((2*L*L) + 1) / (3 * L * L * L));}

static inline int32_t wrapidx(int32_t idx, int32_t n){
        while (idx >= n) {idx -=n;}
return idx;}

// Extract key from the 64-bit value (16-bit key at the lower 16 bits)
static inline uint16_t extract_key(uint64_t value) { return (uint16_t) value;} // Assuming lowest 16 bits only on little endian

static inline void csort64_10(uint64_t** array, size_t n, uint64_t** aux_buffer, size_t* indices) {
    // Clear the indices array (initialized to 0)
    memset(indices, 0, 1024 * sizeof(size_t));

    // Compute the bin sizes (counting phase)
    for (size_t i = 0; i < n; i++) {size_t index_tmp = extract_key((*array)[i]); indices[index_tmp] += 1;}

    // Compute the first index of each bin (prefix sum phase)
    for (size_t i = 1; i < 1024; i++) {indices[i] += indices[i - 1];}

    // Rearrange elements into the aux_buffer based on indices (placement phase)
    for (size_t i = 0; i < n; i++) {
        uint16_t index_tmp = extract_key((*array)[i]);
        size_t pos = --indices[index_tmp];  // Decrease the position before placing the element
        (*aux_buffer)[pos] = (*array)[i];
    }
}

void convolve(kvpair* in, double* temp, double* out, int r, int n) {
    for (int j = 0; j < n; j++) {out[j] = in[j].parts.val;}
    double norm = 1.0 / ((2*r + 1) * (2*r + 1) * (2*r + 1) * (2*r + 1));

    for (int i = 0; i <= 3; i++) {
        float val = 0.0f;
        int idx_hi = r;
        int idx_lo = n - r;
        if (i % 2 == 0) {
            for (int j = n - (r); j < n + (r); j++) {val += out[wrapidx(j, n)];}
            for (int j = 0; j < n; j++) {
                if (idx_hi >= n) {idx_hi = wrapidx(idx_hi, n);}
                if (idx_lo >= n) {idx_lo = wrapidx(idx_lo, n);}
                val += out[idx_hi];
                temp[j] = val;
                val -= out[idx_lo];
                idx_hi += 1; idx_lo += 1;
            }
        } else {
            for (int j = n - (r); j <= n + (r); j++) {val += temp[wrapidx(j, n)];}
            for (int j = 0; j < n; j++) {
                if (idx_hi >= n) {idx_hi = wrapidx(idx_hi, n);}
                if (idx_lo >= n) {idx_lo = wrapidx(idx_lo, n);}
                val += temp[idx_hi];
                out[j] = val;
                val -= temp[idx_lo];
                idx_hi += 1; idx_lo += 1;
            }
        }
    }
    for (int j = 0; j < n; j++) {out[j] *= norm;} // Normalize
}

static inline float get_r(buffer_t *buffer, double freq, float* amp, bool prewhiten){
    double freq_tmp = 1024.0 * freq; //10-bit key
    double min = 0; double max = 0;
    double ref = 0; double res = 0;

    kvpair* input = (kvpair*)(buffer->buf[0]);
    kvpair* sorted = (kvpair*)(buffer->buf[1]);
    double* tmp = (double*)(buffer->buf[0]);
    double* output = (double*)(buffer->buf[2]);
    uint64_t** sort_in = (uint64_t**)(&buffer->buf[0]);
    uint64_t** sort_out = (uint64_t**)(&buffer->buf[1]);

    for (int i = 0; i < buffer->n; i++){
        input[i].parts.key = ((uint16_t)(buffer->x[i] * freq_tmp)) & 0b0000001111111111; // get rid of the 6 most significant bits
        input[i].parts.val = buffer->y[i];
        if (prewhiten){input[i].parts.idx = i;}
    }

    csort64_10(sort_in, buffer->n, sort_out, buffer->pidx); //sort the pairs by phase
    int r = (int)(ceil(0.5 * sqrt((2.0 * (double)(buffer->n)) + 1.0))); r = r >> 1; if (r==0){r=1;} //for 4 convolutions

    convolve(sorted, tmp, output, r, buffer->n);

    double multiplier = 1.0 / (1.0 - corr(r));

    for (int i = 0; i < buffer->n; i++){
        if (amp){
            if (output[i] < min){min = output[i];}
            if (output[i] > max){max = output[i];}
        }
        if (prewhiten) {
            buffer->y[sorted[i].parts.idx] -= output[i];// return 0;
        }
        output[i] = (output[i] - (corr(r) * (double)(sorted[i].parts.val))) * multiplier;
        res += ((double)(sorted[i].parts.val) - output[i]) * ((double)(sorted[i].parts.val) - output[i]);
        ref += ((double)(sorted[i].parts.val) + output[i]) * ((double)(sorted[i].parts.val) + output[i]);
    }

    if (amp){*amp = max - min;}

return (ref / res);}


static inline void binsearch_peak(peak_t *peak, buffer_t *buffer, double df){
    double step = df; float R; double freq_tmp;
    for(int i = 0; i < 12; i++){
        step *= 0.5;

        freq_tmp = peak->freq - step;
        R = get_r(buffer, freq_tmp, NULL, false);
        if(R > peak->r){peak->r = R; peak->freq = freq_tmp;}

        freq_tmp = peak->freq + step;
        R = get_r(buffer, freq_tmp, NULL, false);
        if(R > peak->r){peak->r = R; peak->freq = freq_tmp;}
    }
}

static inline double correct_ihs_res(const double sum, const int n){
    long double logSum = 1.0;
    long double tmp = 1.0; long double denum = 1.0;
    for(int i = 1; i < n; i++){
        denum *= (long double)(i);
        tmp *= sum;
        logSum += tmp / denum;
    }
    return sum - logl(logSum); //asymptotic formula for possible overflow cases may be useful
}

static inline void sortPeaks(peak_t *peaks, int length, buffer_t* buf, int mode, double df, int n) {
    int i, j;
    peak_t key;

    if (mode == 0){
        for (i = 0; i < length; i++){peaks[i].p = correct_ihs_res(peaks[i].p, n);} //Erlang's logarithmic correction'
    }else {
        for (i = 0; i < length; i++){
            if(mode > 1 && mode < 4){
                binsearch_peak(&peaks[i], buf, df);
            }
            if (mode < 4){peaks[i].r = get_r(buf, peaks[i].freq, &peaks[i].amp, false);}
            peaks[i].p = get_z(peaks[i].r, buf->n);
        };

        for (i = 1; i < length; i++) {
            key = peaks[i];
            j = i - 1;
            while (j >= 0 && peaks[j].r < key.r) {
                peaks[j + 1] = peaks[j];
                j = j - 1;
            }
            peaks[j + 1] = key;
        }
    }
};

#define FTEST_H
#endif
