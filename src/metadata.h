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
#include <unistd.h>
#include <errno.h>

#include <fftw3.h>
#include <klib/kvec.h>

#include "utils/common.h"

static inline uint32_t intmax(int32_t a, int32_t b) {return (a > b) ? a : b;}
uint32_t bitCeil(uint32_t n) {
    int exp;
    double mantissa = frexp((double)n, &exp);

    // If the number is a power of two, return it as is
    if (mantissa == 0.5) {return n;}
    // Otherwise, return 2^exp
    return 1 << exp;
}

// Function to check if the path is a file or directory
static int is_directory(const char *path) {
    struct stat path_stat;
    if (stat(path, &path_stat) != 0) {
        perror("stat");
        exit(EXIT_FAILURE);
    }
    return S_ISDIR(path_stat.st_mode);
}

void generate_plans(char *argv[]) {
    // Check if the program is run with root privileges
    if (geteuid() != 0) {
        fprintf(stderr, "Error: %s must be run as root.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    struct stat st;
    if (stat("/opt/ihsnpeaks", &st) == -1) {
        // Directory does not exist, attempt to create it
        if (mkdir("/opt/ihsnpeaks", 0755) == -1) {
            fprintf(stderr, "Error: Failed to create directory '%s': %s\n", "/opt/ihsnpeaks", strerror(errno));
            exit(EXIT_FAILURE);
        }
    }
    printf("Generating FFTW plans - it may take a few minuts.\n");
    fftwf_complex* buffer = fftwf_alloc_complex(1 << 22);
    fftwf_plan plan;
    for (int i = 11; i < 23; i++){
        printf("\rGenerating plan %i out of 12", i - 10); fflush(stdout);
        if (i < FFTW_MEASURE_THRESHOLD) {plan = fftwf_plan_dft_1d(1<<i, buffer, buffer, FFTW_FORWARD, FFTW_PATIENT);}
        else {plan = fftwf_plan_dft_1d(1<<i, buffer, buffer, FFTW_FORWARD, FFTW_MEASURE);}
    }
    int status = fftwf_export_wisdom_to_filename("/opt/ihsnpeaks/plans");
    if (status == 1){printf("Plans successfully saved to /opt/ihsnpeaks/plans");}
    else {printf("Error: Failed to save generated plans");}
    fftwf_free(buffer);
    fftwf_destroy_plan(plan);
    printf("\n");
}

// Function to process the path and populate the target vector
static bool process_path(parameters* params) {
    kv_init(params->targets);
    params->maxLen = 0;

    if (!is_directory(params->target)) {

        // Allocate and populate a target struct with the file path and params pointer
        target_t target;
        target.path = strdup(params->target);
        kv_push(target_t, params->targets, target); // Store the target in the vector

        // Use stat to get the file size
        struct stat file_stat;
        if (stat(params->target, &file_stat) != 0) {
            perror("stat");
        exit(EXIT_FAILURE);
        }

        params->maxSize = file_stat.st_size + 1;

        // Check newline count for this single file
        FILE *file = fopen(params->target, "r");
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

        params->maxSize = file_stat.st_size + 1;
        params->maxLen = newline_count;
        params->gridLen = intmax(1<<11, bitCeil(newline_count * params->gridRatio)); //to be modified (?)
        params->plan = fftwf_plan_dft_1d(params->gridLen, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);
        return true;
    } else {
        DIR *dir = opendir(params->target);
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
                snprintf(full_path, sizeof(full_path), "%s/%s", params->target, entry->d_name);

                // Get the file size
                struct stat file_stat;
                if (stat(full_path, &file_stat) != 0) {
                    perror("stat");
                    exit(EXIT_FAILURE);
                }

                params->avgLen += file_stat.st_size;

                // Check if it's the largest file so far
                if (file_stat.st_size > largest_size) {
                    largest_size = file_stat.st_size;
                    params->maxSize = largest_size + 1;

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

                    params->maxLen = newline_count;
                }

                // Allocate and populate a target struct with the file path and params pointer
                target_t target;
                target.path = strdup(full_path);
                kv_push(target_t, params->targets, target); // Store the target in the vector
            }
        }
        closedir(dir);
        //approximate an near-optimal transform size
        params->avgLen /= DEFAULT_MEASUREMENT_SIZE * kv_size(params->targets);
        params->gridLen = intmax(1<<11, bitCeil((params->avgLen) * (uint64_t)(params->gridRatio)));
        params->plan = fftwf_plan_dft_1d(params->gridLen, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);
    }
return false;}

// Function to free all memory allocated in the kvec and destroy it
static void free_targets(kvec_target_t *targets) {
    for (size_t i = 0; i < kv_size(*targets); i++) {free(kv_A(*targets, i).path);}
    kv_destroy(*targets);
}

#endif // METADATA_H
