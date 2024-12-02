#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <stdint.h>
#include "../include/bsort_key10.h"


//#define __AVX__
#include "../include/simd.h"

//int wrapidx(int idx,int n){return idx % n;}
int32_t wrapidx(int32_t idx, int32_t n){
        while (idx >= n) {idx -=n;}
return idx;}

struct FFT {

    int n;  //number of data points on which convolution is transformed
    int r;
    int half_r;
    __m256 norm;
    __m256* temp;

    // Member function to initialize the FFT struct
    int init(int size) {
        n = size;
        r = int(std::ceil(0.5 * std::sqrt((2.0 * n) + 1))); //for 4 convolutions
        if (r % 2 == 0){r+= 1;}

        float norm_tmp = 1.0;
        for (int i = 0; i <= 3; i++) {norm_tmp /= float(r);}
        norm = _mm256_set1_ps(norm_tmp);
        half_r = (r >> 1);
        temp = (__m256 *) aligned_alloc(64, n * sizeof(__m256));

    return half_r;}




    void convolve(__m256* in, __m256* out) {
        // Copy real part of input to the real part of output (leaving the imaginary part unchanged)
        for (int j = 0; j < n; j++) {out[j] = in[j];}

        for (int i = 0; i <= 3; i++) {
            __m256 val = _mm256_set1_ps(0.0f);
            int idx_hi = half_r;
            int idx_lo = n - half_r;
            if (i % 2 == 0) {
                for (int j = n - (half_r); j < n + (half_r); j++) {val += out[wrapidx(j, n)];}
                for (int j = 0; j < n; j++) {
                    if (idx_hi >= n) {idx_hi = wrapidx(idx_hi, n);}
                    if (idx_lo >= n) {idx_lo = wrapidx(idx_lo, n);}
                    val += out[idx_hi];
                    temp[j] = val;
                    val -= out[idx_lo];
                    idx_hi += 1; idx_lo += 1;
                }
            } else {
                for (int j = n - (half_r); j <= n + (half_r); j++) {val += temp[wrapidx(j, n)];}
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

/*
    void convolve(float* in, float* out) {
        for (int j = 0; j < n; j++) {out[j] = in[j];}

        for (int i = 0; i <= 3; i++){
            float val = 0;
            if (i % 2 == 0) {
                for (int j = n - (half_r); j < n + (half_r); j++) {val += out[wrapidx(j, n)];}
                for (int j = n; j < n << 1; j++) {val += out[wrapidx((j + (half_r)), n)]; temp[wrapidx(j, n)] = val; val -= out[wrapidx((j - (half_r)), n)];}
            }
            else {
                for (int j = n - (half_r); j <= n + (half_r); j++) {val += temp[wrapidx(j, n)];}
                for (int j = n; j < n << 1; j++) {val += temp[wrapidx((j + (half_r)), n)]; out[wrapidx(j, n)] = val; val -= temp[wrapidx((j - (half_r)), n)];}
            }
        }
        for (int j = 0; j < n; j++) {out[j] *= norm;}
    }
*/



        //void convolve(float* y) {convolve(y, y);}

        // Destructor
   void free() {
        ::free(temp);
    }

};
