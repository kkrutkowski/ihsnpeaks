#ifndef SMOOTHER_H
#include "simd.h"
#include "readout.h" //for the "peak_t" struct

static inline void sortPeaks(peak_t *peaks, int length) {
    int i, j;
    peak_t key;

    for (i = 1; i < length; i++) {
        key = peaks[i];
        j = i - 1;
        while (j >= 0 && peaks[j].chi2 > key.chi2) {
            peaks[j + 1] = peaks[j];
            j = j - 1;
        }
        peaks[j + 1] = key;
    }
}

static inline int32_t wrapidx(int32_t idx, int32_t n){
        while (idx >= n) {idx -=n;}
return idx;}














































#define SMOOTHER_H
#endif
