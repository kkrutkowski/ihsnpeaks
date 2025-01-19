#ifndef METADATA_H
#define METADATA_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>

#include "params.h"

#include <fftw3.h>
#include <klib/kvec.h>

static inline uint32_t intmax(int32_t a, int32_t b) {return (a > b) ? a : b;}
uint32_t bitCeil(uint32_t n) {
    int exp;
    double mantissa = frexp((double)n, &exp);

    // If the number is a power of two, return it as is
    if (mantissa == 0.5) {return n;}
    // Otherwise, return 2^exp
    return 1 << exp;
}

// Define the target struct to hold file path and a pointer to the parameters struct
typedef struct {
    char *path;
} target_t;

// Define a kvec type for target_t structs
typedef kvec_t(target_t) kvec_target_t;

// Function to check if the path is a file or directory
static int is_directory(const char *path) {
    struct stat path_stat;
    if (stat(path, &path_stat) != 0) {
        perror("stat");
        exit(EXIT_FAILURE);
    }
    return S_ISDIR(path_stat.st_mode);
}

// Function to process the path and populate the target vector
static bool process_path(char* path, kvec_target_t *targets, uint32_t* maxLen, uint32_t* maxSize, uint64_t* avgLen, uint32_t* gridLen, uint32_t* gridRatio, fftwf_plan* plan) {
    kv_init(*targets);
    *maxLen = 0;

    if (!is_directory(path)) {

        // Allocate and populate a target struct with the file path and params pointer
        target_t target;
        target.path = strdup(path);
        kv_push(target_t, *targets, target); // Store the target in the vector

        // Use stat to get the file size
        struct stat file_stat;
        if (stat(path, &file_stat) != 0) {
            perror("stat");
        exit(EXIT_FAILURE);
        }

        *maxSize = file_stat.st_size + 1;

        // Check newline count for this single file
        FILE *file = fopen(path, "r");
        if (!file) {
            perror("fopen");
            exit(EXIT_FAILURE);
        }

        uint32_t newline_count = 0;
        int ch;
        while ((ch = fgetc(file)) != EOF) {
            if (ch == '\n') {
                newline_count++;
            }
        }
        fclose(file);

        *maxSize = file_stat.st_size + 1;
        *maxLen = newline_count;
        *gridLen = intmax(1<<11, bitCeil(newline_count * *gridRatio)); //to be modified (?)
        *plan = fftwf_plan_dft_1d(*gridLen, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);
        return true;
    } else {
        DIR *dir = opendir(path);
        if (!dir) {
            perror("opendir");
            exit(EXIT_FAILURE);
        }

        struct dirent *entry;
        off_t largest_size = 0;

        while ((entry = readdir(dir)) != NULL) {
            if (entry->d_type == DT_REG) {
                // Construct the full path for each file
                char full_path[256];
                snprintf(full_path, sizeof(full_path), "%s/%s", path, entry->d_name);

                // Get the file size
                struct stat file_stat;
                if (stat(full_path, &file_stat) != 0) {
                    perror("stat");
                    exit(EXIT_FAILURE);
                }

                *avgLen += file_stat.st_size;

                // Check if it's the largest file so far
                if (file_stat.st_size > largest_size) {
                    largest_size = file_stat.st_size;
                    *maxSize = largest_size + 1;

                    // Open and count newlines in this file
                    FILE *file = fopen(full_path, "r");
                    if (!file) {
                        perror("fopen");
                        exit(EXIT_FAILURE);
                    }

                    uint32_t newline_count = 0;
                    int ch;
                    while ((ch = fgetc(file)) != EOF) {
                        if (ch == '\n') {
                            newline_count++;
                        }
                    }
                    fclose(file);

                    *maxLen = newline_count;
                }

                // Allocate and populate a target struct with the file path and params pointer
                target_t target;
                target.path = strdup(full_path);
                kv_push(target_t, *targets, target); // Store the target in the vector
            }
        }
        closedir(dir);
        //approximate an near-optimal transform size
        *avgLen /= DEFAULT_MEASUREMENT_SIZE * kv_size(*targets);
        *gridLen = intmax(1<<11, bitCeil((*avgLen) * (uint64_t)(*gridRatio)));
        *plan = fftwf_plan_dft_1d(*gridLen, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);
    }

return false;}

// Function to free all memory allocated in the kvec and destroy it
static void free_targets(kvec_target_t *targets) {
    for (size_t i = 0; i < kv_size(*targets); i++) {
        free(kv_A(*targets, i).path);
    }
    kv_destroy(*targets);
}

#endif // METADATA_H
