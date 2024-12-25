#ifndef SORT_H
#define SORT_H

#include <stdint.h>

typedef struct {uint32_t freq; float power;} peak;

#define INSERTION_SORT_THRESHOLD 32

// Insertion sort for small partitions
void insertion_sort(int64_t *arr, int64_t left, int64_t right) {
    for (int64_t i = left + 1; i <= right; i++) {
        int64_t key = arr[i];
        int64_t j = i - 1;
        while (j >= left && arr[j] < key) {  // Sort in descending order
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

// Swap two elements in the array
void swap(int64_t *a, int64_t *b) {
    int64_t temp = *a;
    *a = *b;
    *b = temp;
}

// Partition the array using two pivots
int64_t partition(int64_t *arr, int64_t low, int64_t high, int64_t *lp) {
    // Ensure the first and last elements are larger than the pivots
    if (arr[low] < arr[high]) {swap(&arr[low], &arr[high]);}

    int64_t p = arr[low], q = arr[high];  // Two pivots
    int64_t j = low + 1;  // Left boundary for values greater than p
    int64_t g = high - 1; // Right boundary for values smaller than q
    int64_t k = low + 1;  // Current index to scan through the array

    while (k <= g) {
        if (arr[k] > p) {
            swap(&arr[k], &arr[j]);  // Elements greater than p go to the left side
            j++;
        } else if (arr[k] < q) {
            swap(&arr[k], &arr[g]);  // Elements smaller than q go to the right side
            g--;
            if (arr[k] > p) {
                swap(&arr[k], &arr[j]);
                j++;
            }
        }
        k++;
    }

    // Place pivots in correct positions
    j--;
    g++;

    swap(&arr[low], &arr[j]);
    swap(&arr[high], &arr[g]);

    *lp = j;  // Left pivot index
    return g; // Right pivot index
}

// Quicksort with tail recursion and insertion sort for small partitions
void quicksort(int64_t *arr, int64_t low, int64_t high) {
    while (low < high) {
        if (high - low + 1 <= INSERTION_SORT_THRESHOLD) {
            insertion_sort(arr, low, high);
            break;
        }

        int64_t lp, rp;
        rp = partition(arr, low, high, &lp);

        // Recursively process the two subarrays defined by the pivots
        quicksort(arr, low, lp - 1);       // Left partition
        quicksort(arr, rp + 1, high);      // Right partition

        low = lp + 1;                     // Set the boundary for the next iteration
        high = rp - 1;
    }
}

#endif
