#ifndef FTEST_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

#include <fdist.h>

#include "simd.h"

#ifndef BUFFER_T
#define BUFFER_T
typedef struct {
    bool allocated;
    uint8_t loc_iter;
    uint32_t   len;
    uint32_t     n;
    uint32_t     r;
    uint32_t terms;
    uint32_t memBlockSize;

    peak_t* peaks;
    uint32_t nPeaks;

    char*  readBuf;
    double*      x;
    float*       y;
    float*      dy;

    size_t* pidx; //phase indices for counting sort

    void** buf;

    float magnitude;
    uint32_t gridSize; uint32_t nGrids;

    complex float ** grids; //used to compute the FFT
    uint32_t** gidx; //grid indices for NFFT
    float** gdist;
    float** weights;

    sds spectrum;
    sds outBuf;
} buffer_t;
#endif

#ifndef PEAK_T
#define PEAK_T
typedef struct {
    double freq;
    float  p;
    float  amp;
    float  r;
} peak_t;
#endif

#ifndef KVPAIR
#define KVPAIR
typedef union {
    uint64_t data; // The raw 64-bit representation
    struct {
        uint16_t key;        // 16-bit key
        uint16_t idx; // 16-bit index
        float val;         // 32-bit floating-point value
    } parts;
} kvpair;
#endif

static constexpr double corr[128] = {
0.2345679012345679, 0.136, 0.09620991253644315, 0.0745313214449017,
0.06085649887302779, 0.051433773327264454, 0.04454320987654321, 0.039283533482597194,
0.03513631724741216, 0.03178202497930389, 0.029012903756061477, 0.026688,
0.02470829311249979, 0.023002173110828653, 0.021516565405659428, 0.020211295693389357,
0.01905539358600583, 0.01802459874044973, 0.017099636429024987, 0.016264999056891223,
0.015508068471958444, 0.014818472793781436, 0.01418760775550697, 0.013608275463454852,
0.013074408284395394, 0.012580855336956012, 0.012123215627347859, 0.011697706356791744,
0.011301058043909065, 0.010930430300333508, 0.010583343664724364, 0.010257624032771962,
0.009951357048573129, 0.009662850434336146, 0.009390602691730626, 0.009133276951906986,
0.00888967901234568, 0.008658738798728678, 0.008439494644439893, 0.008231079900371505,
0.00803271148172309, 0.00784368003256666, 0.007663341447697782, 0.007491109538149904,
0.0073264496643315625, 0.00716887319104991, 0.007017932643242455, 0.006873217461237486,
0.00673435027072411, 0.00660098359605591, 0.006472796956604898, 0.006349494295072527,
0.006230801694307874, 0.006116465345563689, 0.00600624973646644, 0.005899936031470022,
0.00579732062135284, 0.005698213821524509, 0.005602438701629935, 0.005509830031254922,
0.005420233328514791, 0.005333504, 0.00524950656200525, 0.005168113934218386,
0.005089206798123386, 0.005012673013303795, 0.004938407085640739, 0.004866309682101212,
0.0047962871874230855, 0.004728251298535982, 0.0046621186530228535, 0.004597810488334905,
0.004535252328830401, 0.0044743736980225515, 0.004415107853698761, 0.004357391543818082,
0.004301164781309792, 0.004246370636087354, 0.0041929550427617064, 0.00414086662268848,
0.004090056519117834, 0.004040478244334994, 0.003992087537786129, 0.0039448422342794175,
0.0038987021414362984, 0.0038536289256442557, 0.0038095860058309037, 0.003766538454440658,
0.0037244529050505996, 0.0036832974661119263, 0.0036430416403483165, 0.0036036562493830574,
0.003565113363203436, 0.0035273862341040378, 0.003490449234780609, 0.003454277800273388,
0.003418848373483482, 0.0033841383540083652, 0.003350126050062946, 0.003316790633271273,
0.0032841120961308967, 0.003252071211967325, 0.0032206494972101597, 0.0031898291758353766,
0.0031595931458300044, 0.0031299249475462746, 0.003100808733822179, 0.003072229241754485,
0.0030441717660185403, 0.0030166221336368885, 0.0029895666801057247, 0.0029629922267946957,
0.002936886059541511, 0.0029112359083683062, 0.0028860299282517754, 0.002861256680883749,
0.002836905117363205, 0.0028129645617646906, 0.0027894246955318153, 0.002766275542647873,
0.002743507455538832, 0.0027211111016668225, 0.002699077450775, 0.002677397762747144,
0.0026560635760477084, 0.002635066696710202, 0.0026143991878437908, 0.002594053359629889};


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

static inline float get_r(buffer_t *buffer, double freq, float* amp){
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
    }

    csort64_10(sort_in, buffer->n, sort_out, buffer->pidx); //sort the pairs by phase
    int r = (int)(ceil(0.5 * sqrt((2.0 * (double)(buffer->n)) + 1.0))); r = r >> 1; if (r==0){r=1;} //for 4 convolutions

    convolve(sorted, tmp, output, r, buffer->n);

    double multiplier = 1.0 / (1.0 - corr[r-1]);

    for (int i = 0; i < buffer->n; i++){
        if (amp){
            if (output[i] < min){min = output[i];}
            if (output[i] > max){max = output[i];}
        }
        output[i] = (output[i] - (corr[r-1] * (double)(sorted[i].parts.val))) * multiplier;
        res += ((double)(sorted[i].parts.val) - output[i]) * ((double)(sorted[i].parts.val) - output[i]);
        ref += ((double)(sorted[i].parts.val) + output[i]) * ((double)(sorted[i].parts.val) + output[i]);
    }

    if (amp){*amp = max - min;}

return (ref / res);}


static inline void binsearch_peak(peak_t *peak, buffer_t *buffer, double df){
    double step = df; float R; double freq_tmp;
    for(int i = 0; i < 10; i++){
        step *= 0.5;

        freq_tmp = peak->freq - step;
        R = get_r(buffer, freq_tmp, NULL);
        if(R > peak->r){peak->r = R; peak->freq = freq_tmp;}

        freq_tmp = peak->freq + step;
        R = get_r(buffer, freq_tmp, NULL);
        if(R > peak->r){peak->r = R; peak->freq = freq_tmp;}
    }
}


static inline void sortPeaks(peak_t *peaks, int length, buffer_t* buf, int mode, double df) {
    int i, j;
    peak_t key;

    for (i = 0; i < length; i++){
        if(mode > 1 && mode < 4){binsearch_peak(&peaks[i], buf, df);}
        peaks[i].r = get_r(buf, peaks[i].freq, &peaks[i].amp);
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

#define FTEST_H
#endif
