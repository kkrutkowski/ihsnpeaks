#ifndef SMOOTHER_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

#include "simd.h"
#include "readout.h" //for the "peak_t" struct

typedef uint64_t T; // 64-bit key-value pair
typedef size_t U;


#ifndef KVPAIR
#define KVPAIR
typedef union {
    uint64_t data; // The raw 64-bit representation
    struct {
        uint16_t key;        // 16-bit key
        uint16_t placeholder; // 16-bit placeholder
        float value;         // 32-bit floating-point value
    } parts;
} kvpair;
#endif


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

// Extract key from the 64-bit value (16-bit key at the lower 16 bits)
static inline uint16_t extract_key(uint64_t value) { return (uint16_t) value; } // Assuming lowest 16 bits only on little endian

static inline void bsort64_10(uint64_t** array, size_t n, uint64_t* aux_buffer, size_t* indices) {
    // Clear the indices array (initialized to 0)
    memset(indices, 0, 1024 * sizeof(size_t));

    // Compute the bin sizes (counting phase)
    for (size_t i = 0; i < n; i++) {size_t index_tmp = extract_key((*array)[i]); indices[index_tmp] += 1;}

    // Compute the first index of each bin (prefix sum phase)
    //to be vectorized - https://stackoverflow.com/questions/69694914/accumulating-a-running-total-prefix-sum-horizontally-across-an-m256i-vector
    for (size_t i = 1; i < 1024; i++) {indices[i] += indices[i - 1];}

    // Rearrange elements into the aux_buffer based on indices (placement phase)
    for (size_t i = 0; i < n; i++) {
        uint16_t index_tmp = extract_key((*array)[i]);
        size_t pos = --indices[index_tmp];  // Decrease the position before placing the element
        (aux_buffer)[pos] = (*array)[i];
    }

    // Swap the pointers
    uint64_t* tmp = aux_buffer;
    aux_buffer = *array;
    *array = tmp;
}












































#define SMOOTHER_H
#endif
