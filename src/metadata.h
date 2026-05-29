#ifndef METADATA_H
#define METADATA_H

#include <dirent.h>
#include <klib/kvec.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "scaling.h"
#include "utils/common.h"

static inline uint32_t intmax(int32_t a, int32_t b) { return (a > b) ? a : b; }
uint32_t bitCeil(uint32_t n) {
    int exp;
    double mantissa = frexp((double)n, &exp);

    // If the number is a power of two, return it as is
    if (mantissa == 0.5) {
        return n;
    }
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
    (void)argv;
    printf("FFTW runtime plan generation is obsolete: ihsnpeaks now uses the built-in PSWF NuFFT backend.\n");
}

static inline size_t round_metadata_alloc(size_t size) { return (size + 63U) & ~(size_t)63U; }

static void inspect_dat_file(const char *path, uint32_t *line_count, double *time_span) {
    FILE *file = fopen(path, "r");
    if (!file) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    double x = 0.0, y = 0.0, dy = 0.0;
    double first_x = 0.0, last_x = 0.0;
    uint32_t count = 0;
    while (fscanf(file, "%lf %lf %lf", &x, &y, &dy) == 3) {
        if (count == 0) first_x = x;
        last_x = x;
        ++count;
    }
    fclose(file);

    *line_count = count;
    *time_span = count > 1 ? last_x - first_x : 0.0;
}

static uint32_t estimate_frequency_count(const parameters *params, double time_span) {
    if (time_span <= 0.0 || params->fmax <= params->fmin) return 1;
    double count;
    if (mode_uses_direct_gb_grid(params->mode)) {
        double fstep = gb_direct_frequency_step(params->maxLen, time_span, params->oversamplingFactor, params->gbAlpha);
        if (fstep <= 0.0) return 1;
        count = ((double)params->fmax - (double)params->fmin) / fstep;
    } else {
        count = ((double)params->fmax - (double)params->fmin) * (double)params->nterms * (double)params->oversamplingFactor * time_span * 0.5;
    }
    if (count > (double)UINT32_MAX - 2.0) return UINT32_MAX - 1U;
    return (uint32_t)floor(count) + 1U;
}

static int pswf_backend(nufft1_mode mode) { return mode == NUFFT1_PSWF21 ? 2 : 1; }

static void initialize_nufft_plan(parameters *params) {
    double beta = params->gridMode == NUFFT1_PSWF21 ? F_PSWF21_BETA : F_PSWF43_BETA;
    double gamma = params->gridMode == NUFFT1_PSWF21 ? F_PSWF21_GAMMA : F_PSWF43_GAMMA;
    int base_len = nufft1_optimize_plan_size((int)params->maxFreqCount, (int)params->maxLen, 0, F_ALPHA, beta, gamma, pswf_backend(params->gridMode));
    if (params->gridMode == NUFFT1_PSWF43) {
        base_len = nufft1_pswf43_plan_len_from_base(base_len);
    }
    if (base_len < (1 << 8)) base_len = 1 << 8;
    if (base_len > (1 << 20)) base_len = 1 << 20;
    if (params->gridMode == NUFFT1_PSWF43) base_len = nufft1_pswf43_plan_len_from_base(base_len);

    params->gridLen = (uint32_t)base_len;
    params->outputLen = params->gridMode == NUFFT1_PSWF43 ? (uint32_t)nufft1_pswf43_output_len_for_plan(base_len) : (uint32_t)base_len;
    params->ladderLevels = (uint32_t)nufft1_twiddle_ladder_levels((int)params->maxFreqCount, (int)params->outputLen);
    if (params->ladderLevels > 8U) params->ladderLevels = 8U;
    if (params->ladderLevels == 0U) params->ladderLevels = 1U;

    params->nufftExternalSizes = nufft1_external_sizes_for_plan((int)params->gridLen, params->gridMode);
    params->nufftTwiddleReal = aligned_alloc(64, round_metadata_alloc(params->nufftExternalSizes.twiddle_len * sizeof(float)));
    params->nufftTwiddleImag = aligned_alloc(64, round_metadata_alloc(params->nufftExternalSizes.twiddle_len * sizeof(float)));
    if (!params->nufftTwiddleReal || !params->nufftTwiddleImag) {
        fprintf(stderr, "Failed to allocate shared NuFFT twiddles\n");
        exit(EXIT_FAILURE);
    }
    params->nufftPlan = nufft1_initialize_shared((int)params->maxLen, (int)params->gridLen, params->nterms, params->gridMode, params->nufftTwiddleReal,
                                                 params->nufftTwiddleImag);
    if (!params->nufftPlan) {
        fprintf(stderr, "Failed to initialize shared NuFFT plan\n");
        exit(EXIT_FAILURE);
    }
}

// Function to process the path and populate the target vector
static bool process_path(parameters *params) {
    kv_init(params->targets);
    params->maxLen = 0;
    params->maxTimeSpan = 0.0;

    if (!is_directory(params->target)) {
        // Allocate and populate a target struct with the file path and params pointer
        target_t target;
        target.path = strdup(params->target);
        kv_push(target_t, params->targets, target);  // Store the target in the vector

        // Use stat to get the file size
        struct stat file_stat;
        if (stat(params->target, &file_stat) != 0) {
            perror("stat");
            exit(EXIT_FAILURE);
        }

        params->maxSize = file_stat.st_size + 1;

        params->maxSize = file_stat.st_size + 1;
        uint32_t newline_count = 0;
        double time_span = 0.0;
        inspect_dat_file(params->target, &newline_count, &time_span);
        params->maxLen = newline_count;
        params->maxTimeSpan = time_span;
        params->maxFreqCount = estimate_frequency_count(params, params->maxTimeSpan);
        initialize_nufft_plan(params);
        return true;
    } else {
        DIR *dir = opendir(params->target);
        if (!dir) {
            perror("opendir");
            exit(EXIT_FAILURE);
        }

        struct dirent *entry;
        off_t largest_size = 0;
        char largest_path[256] = {0};

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
                    snprintf(largest_path, sizeof(largest_path), "%s", full_path);
                }

                // Allocate and populate a target struct with the file path and params pointer
                target_t target;
                target.path = strdup(full_path);
                kv_push(target_t, params->targets, target);  // Store the target in the vector
            }
        }
        closedir(dir);

        if (largest_path[0] == '\0') {
            fprintf(stderr, "No regular files found in '%s'\n", params->target);
            exit(EXIT_FAILURE);
        }

        uint32_t newline_count = 0;
        double time_span = 0.0;
        inspect_dat_file(largest_path, &newline_count, &time_span);
        params->maxLen = newline_count;
        params->maxTimeSpan = time_span;

        // approximate an near-optimal transform size
        params->avgLen /= DEFAULT_MEASUREMENT_SIZE * kv_size(params->targets);
        params->maxFreqCount = estimate_frequency_count(params, params->maxTimeSpan);
        initialize_nufft_plan(params);
    }
    return false;
}

// Function to free all memory allocated in the kvec and destroy it
static void free_targets(kvec_target_t *targets) {
    for (size_t i = 0; i < kv_size(*targets); i++) {
        free(kv_A(*targets, i).path);
    }
    kv_destroy(*targets);
}

#endif  // METADATA_H
