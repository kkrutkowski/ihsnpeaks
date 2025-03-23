#ifndef FTEST_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

#include <fdist.h>

#include "common.h"
#include "simd.h"

// Expected values of self-correlation of the convolved array elements
static constexpr double corr[128] = {
0.23456790, 0.13600000, 0.09620991, 0.07453132, 0.06085650, 0.05143377, 0.04454321, 0.03928353, 
0.03513632, 0.03178202, 0.02901290, 0.02668800, 0.02470829, 0.02300217, 0.02151657, 0.02021130, 
0.01905539, 0.01802460, 0.01709964, 0.01626500, 0.01550807, 0.01481847, 0.01418761, 0.01360828, 
0.01307441, 0.01258086, 0.01212322, 0.01169771, 0.01130106, 0.01093043, 0.01058334, 0.01025762, 
0.00995136, 0.00966285, 0.00939060, 0.00913328, 0.00888968, 0.00865874, 0.00843949, 0.00823108, 
0.00803271, 0.00784368, 0.00766334, 0.00749111, 0.00732645, 0.00716887, 0.00701793, 0.00687322, 
0.00673435, 0.00660098, 0.00647280, 0.00634949, 0.00623080, 0.00611647, 0.00600625, 0.00589994, 
0.00579732, 0.00569821, 0.00560244, 0.00550983, 0.00542023, 0.00533350, 0.00524951, 0.00516811, 
0.00508921, 0.00501267, 0.00493841, 0.00486631, 0.00479629, 0.00472825, 0.00466212, 0.00459781, 
0.00453525, 0.00447437, 0.00441511, 0.00435739, 0.00430116, 0.00424637, 0.00419296, 0.00414087, 
0.00409006, 0.00404048, 0.00399209, 0.00394484, 0.00389870, 0.00385363, 0.00380959, 0.00376654, 
0.00372445, 0.00368330, 0.00364304, 0.00360366, 0.00356511, 0.00352739, 0.00349045, 0.00345428, 
0.00341885, 0.00338414, 0.00335013, 0.00331679, 0.00328411, 0.00325207, 0.00322065, 0.00318983, 
0.00315959, 0.00312992, 0.00310081, 0.00307223, 0.00304417, 0.00301662, 0.00298957, 0.00296299, 
0.00293689, 0.00291124, 0.00288603, 0.00286126, 0.00283691, 0.00281296, 0.00278942, 0.00276628, 
0.00274351, 0.00272111, 0.00269908, 0.00267740, 0.00265606, 0.00263507, 0.00261440, 0.00259405
};
// Expected original to smoothened variance ratios, post-decorellation
static constexpr double varr[128] = {
5.15281501,   9.42469295,   13.64084054,  17.83787926,  22.02605065,  26.20938269,  30.38978405,  34.56827624, 
38.74545543,  42.92169298,  47.09723226,  51.27223944,  55.44683184,  59.62109465,  63.79509122,  67.96886966, 
72.14246713,  76.31591282,  80.48922993,  84.66243721,  88.83554989,  93.00858053,  97.18153956,  101.35443571, 
105.52727634, 109.70006772, 113.87281520, 118.04552338, 122.21819624, 126.39083725, 130.56344943, 134.73603544, 
138.90859760, 143.08113799, 147.25365844, 151.42616059, 155.59864589, 159.77111565, 163.94357106, 168.11601317, 
172.28844294, 176.46086124, 180.63326886, 184.80566651, 188.97805486, 193.15043450, 197.32280598, 201.49516981, 
205.66752644, 209.83987631, 214.01221980, 218.18455729, 222.35688910, 226.52921555, 230.70153692, 234.87385350, 
239.04616552, 243.21847322, 247.39077681, 251.56307651, 255.73537250, 259.90766496, 264.07995405, 268.25223993, 
272.42452275, 276.59680264, 280.76907974, 284.94135417, 289.11362604, 293.28589546, 297.45816254, 301.63042737, 
305.80269004, 309.97495063, 314.14720924, 318.31946594, 322.49172081, 326.66397390, 330.83622530, 335.00847506, 
339.18072324, 343.35296991, 347.52521510, 351.69745889, 355.86970131, 360.04194242, 364.21418225, 368.38642086, 
372.55865828, 376.73089455, 380.90312972, 385.07536381, 389.24759685, 393.41982890, 397.59205996, 401.76429009, 
405.93651929, 410.10874761, 414.28097507, 418.45320169, 422.62542750, 426.79765252, 430.96987678, 435.14210030, 
439.31432310, 443.48654519, 447.65876661, 451.83098736, 456.00320746, 460.17542694, 464.34764581, 468.51986409, 
472.69208178, 476.86429892, 481.03651550, 485.20873155, 489.38094708, 493.55316211, 497.72537663, 501.89759068, 
506.06980425, 510.24201737, 514.41423004, 518.58644227, 522.75865407, 526.93086546, 531.10307644, 535.27528702
};



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
    for(int i = 0; i < 12; i++){
        step *= 0.5;

        freq_tmp = peak->freq - step;
        R = get_r(buffer, freq_tmp, NULL);
        if(R > peak->r){peak->r = R; peak->freq = freq_tmp;}

        freq_tmp = peak->freq + step;
        R = get_r(buffer, freq_tmp, NULL);
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
            if(mode > 1 && mode < 4){binsearch_peak(&peaks[i], buf, df);}
            if (mode < 4){peaks[i].r = get_r(buf, peaks[i].freq, &peaks[i].amp);}
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
