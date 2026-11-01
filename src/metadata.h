#ifndef METADATA_H
#define METADATA_H

#include <dirent.h>
#include <errno.h>
#include <klib/kvec.h>
#include <math.h>
#include <qfits/qfits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "scaling.h"
#include "utils/common.h"
#include "utils/convolution.h"

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
        const eval_method_t *method = eval_method_for_params(params);
        double fstep = method->direct_frequency_step(params->maxLen, time_span, params);
        if (fstep <= 0.0) return 1;
        count = ((double)params->fmax - (double)params->fmin) / fstep;
    } else {
        count = ((double)params->fmax - (double)params->fmin) * (double)params->nterms * (double)params->oversamplingFactor * time_span * 0.5;
    }
    if (count > (double)UINT32_MAX - 2.0) return UINT32_MAX - 1U;
    return (uint32_t)floor(count) + 1U;
}

static int pswf_backend(nufft1_mode mode) { return mode == NUFFT1_PSWF21 ? 2 : 1; }

static bool has_dat_suffix(const char *path) {
    size_t len = strlen(path);
    return len >= 4U && strcmp(path + len - 4U, ".dat") == 0;
}

static bool has_fits_suffix(const char *path) {
    size_t len = strlen(path);
    return (len >= 5U && strcmp(path + len - 5U, ".fits") == 0) || (len >= 4U && strcmp(path + len - 4U, ".fit") == 0);
}

static void inspect_fits_file(const char *path, uint32_t *line_count, double *time_span) {
    qfits_table *table = qfits_table_open(path, 1);
    if (!table) {
        *line_count = 0;
        *time_span = 0.0;
        return;
    }

    int col_time = -1;
    for (int i = 0; i < table->nc; i++) {
        if (strcmp(table->col[i].tlabel, "TIME") == 0) {
            col_time = i;
            break;
        }
    }

    *line_count = (uint32_t)table->nr;

    if (col_time >= 0) {
        double null_time = 0.0;
        void *time_data = qfits_query_column_data(table, col_time, NULL, &null_time);
        if (time_data) {
            tfits_type t_type = table->col[col_time].atom_type;
            double first_x = 0.0, last_x = 0.0;
            int valid_count = 0;
            for (int i = 0; i < table->nr; i++) {
                double t = 0.0;
                if (t_type == TFITS_BIN_TYPE_D)
                    t = ((double *)time_data)[i];
                else if (t_type == TFITS_BIN_TYPE_E)
                    t = ((float *)time_data)[i];

                if (!isnan(t) && !isinf(t) && t != null_time) {
                    if (valid_count == 0) first_x = t;
                    last_x = t;
                    valid_count++;
                }
            }
            *time_span = valid_count > 1 ? last_x - first_x : 0.0;
            free(time_data);
        } else {
            *time_span = 0.0;
        }
    } else {
        *time_span = 0.0;
    }

    qfits_table_close(table);
}

static uint32_t nufft_output_len_for_grid(uint32_t grid_len, nufft1_mode mode) {
    return mode == NUFFT1_PSWF43 ? (uint32_t)nufft1_pswf43_output_len_for_plan((int)grid_len) : grid_len;
}

static void free_nufft_plan_cache(parameters *params) {
    if (!params || !params->nufftPlanCache) return;
    for (uint32_t i = 0; i < params->nufftPlanCount; ++i) {
        nufft1_free_plan(params->nufftPlanCache[i].plan);
    }
    free(params->nufftPlanCache);
    free(params->nufftTwiddleReal);
    free(params->nufftTwiddleImag);
    params->nufftPlanCache = NULL;
    params->nufftTwiddleReal = NULL;
    params->nufftTwiddleImag = NULL;
    params->nufftPlanCount = 0;
    params->nufftExternalSizes = (nufft1_external_sizes){0};
}

static void initialize_nufft_plan(parameters *params) {
    double beta = params->gridMode == NUFFT1_PSWF21 ? F_PSWF21_BETA : F_PSWF43_BETA;
    double gamma = params->gridMode == NUFFT1_PSWF21 ? F_PSWF21_GAMMA : F_PSWF43_GAMMA;
    int plan_degree = periodogram_uses_aov(params->periodogramMethod) ? params->nterms : 0;
    int base_len = nufft1_optimize_plan_size((int)params->maxFreqCount, (int)params->maxLen, plan_degree, F_ALPHA, beta, gamma, pswf_backend(params->gridMode));
    if (params->gridMode == NUFFT1_PSWF43) {
        base_len = nufft1_pswf43_plan_len_from_base(base_len);
    }
    if (base_len < (1 << 8)) base_len = 1 << 8;
    if (base_len > (1 << 20)) base_len = 1 << 20;
    if (params->gridMode == NUFFT1_PSWF43) base_len = nufft1_pswf43_plan_len_from_base(base_len);

    params->gridLen = (uint32_t)base_len;
    params->outputLen = nufft_output_len_for_grid(params->gridLen, params->gridMode);
    params->nufftExternalSizes = nufft1_external_sizes_for_plan((int)params->gridLen, params->gridMode);
    int nufft_factors = periodogram_uses_aov(params->periodogramMethod) ? 2 * params->nterms : params->nterms;

    uint32_t plan_count = 0;
    for (uint32_t len = 1U << 8; len < params->gridLen; len <<= 1U) ++plan_count;
    ++plan_count;

    params->nufftPlanCache = calloc(plan_count, sizeof(*params->nufftPlanCache));
    if (!params->nufftPlanCache) {
        fprintf(stderr, "Failed to allocate shared NuFFT plan cache\n");
        exit(EXIT_FAILURE);
    }

    uint32_t idx = 0;
    for (uint32_t len = 1U << 8; len < params->gridLen; len <<= 1U) {
        params->nufftPlanCache[idx++].gridLen = len;
    }
    params->nufftPlanCache[idx++].gridLen = params->gridLen;
    params->nufftPlanCount = idx;

    /* The fused-Radix-2 twiddle layout is not a suffix of the next-larger
     * plan's layout, so every cached plan needs its own region inside the
     * shared allocation.  The total is dominated by the largest plan and is
     * roughly twice the old suffix-shared size. */
    size_t total_twiddle_len = 0;
    for (uint32_t i = 0; i < params->nufftPlanCount; ++i) {
        nufft_plan_cache_entry_t *entry = &params->nufftPlanCache[i];
        if (params->gridMode == NUFFT1_PSWF43) entry->gridLen = (uint32_t)nufft1_pswf43_plan_len_from_base((int)entry->gridLen);
        entry->externalSizes = nufft1_external_sizes_for_plan((int)entry->gridLen, params->gridMode);
        total_twiddle_len += entry->externalSizes.twiddle_len;
    }
    params->nufftTwiddleReal = aligned_alloc(64, round_metadata_alloc(total_twiddle_len * sizeof(float)));
    params->nufftTwiddleImag = aligned_alloc(64, round_metadata_alloc(total_twiddle_len * sizeof(float)));
    if (!params->nufftTwiddleReal || !params->nufftTwiddleImag) {
        fprintf(stderr, "Failed to allocate shared NuFFT twiddles\n");
        free_nufft_plan_cache(params);
        exit(EXIT_FAILURE);
    }

    size_t twiddle_offset = 0;
    for (uint32_t i = 0; i < params->nufftPlanCount; ++i) {
        nufft_plan_cache_entry_t *entry = &params->nufftPlanCache[i];
        entry->outputLen = nufft_output_len_for_grid(entry->gridLen, params->gridMode);
        entry->twiddleReal = params->nufftTwiddleReal + twiddle_offset;
        entry->twiddleImag = params->nufftTwiddleImag + twiddle_offset;
        twiddle_offset += entry->externalSizes.twiddle_len;
        /* nufft1_initialize_shared generates the fused layout into this
         * plan's own region; the old suffix-sharing path is invalid. */
        entry->plan =
            nufft1_initialize_shared((int)params->maxLen, (int)entry->gridLen, nufft_factors, params->gridMode, entry->twiddleReal, entry->twiddleImag);
        if (!entry->plan) {
            fprintf(stderr, "Failed to initialize shared NuFFT plan\n");
            free_nufft_plan_cache(params);
            exit(EXIT_FAILURE);
        }
    }
}

/* ---------- Input target discovery (file / directory / .tsv list) ---------- */

static bool has_tsv_suffix(const char *path) {
    size_t len = strlen(path);
    return len >= 4U && strcmp(path + len - 4U, ".tsv") == 0;
}

static bool has_photometry_suffix(const char *path) { return has_dat_suffix(path) || has_fits_suffix(path); }

/* Prepend the optional --prefix to a photometry path. Absolute paths bypass the
 * prefix (there is no sensible way to prefix them); *bypassed reports that. */
static char *resolve_photometry_path(const parameters *params, const char *path, bool *bypassed) {
    if (bypassed) *bypassed = false;
    if (!params->prefix || params->prefix[0] == '\0' || path[0] == '\0') return strdup(path);
    if (path[0] == '/') {
        if (bypassed) *bypassed = true;
        return strdup(path);
    }

    size_t prefix_len = strlen(params->prefix);
    size_t path_len = strlen(path);
    bool need_slash = params->prefix[prefix_len - 1U] != '/';
    char *resolved = malloc(prefix_len + path_len + (need_slash ? 2U : 1U));
    if (!resolved) return NULL;
    memcpy(resolved, params->prefix, prefix_len);
    size_t offset = prefix_len;
    if (need_slash) resolved[offset++] = '/';
    memcpy(resolved + offset, path, path_len + 1U);
    return resolved;
}

static char *join_directory_entry(const char *dir, const char *name) {
    size_t dir_len = strlen(dir);
    size_t name_len = strlen(name);
    char *joined = malloc(dir_len + name_len + 2U);
    if (!joined) return NULL;
    memcpy(joined, dir, dir_len);
    size_t offset = dir_len;
    joined[offset++] = '/';
    memcpy(joined + offset, name, name_len + 1U);
    return joined;
}

/* Batch sizing accumulators. Only the largest file per format is inspected. */
typedef struct {
    off_t largest_size;
    off_t largest_dat_size;
    off_t largest_fits_size;
    size_t file_count;
    uint64_t total_bytes;
    char *largest_dat_path;
    char *largest_fits_path;
} batch_meta_t;

static void batch_meta_init(batch_meta_t *meta) { memset(meta, 0, sizeof(*meta)); }

static void batch_meta_add(batch_meta_t *meta, const char *path, const struct stat *file_stat) {
    meta->file_count++;
    meta->total_bytes += (uint64_t)file_stat->st_size;
    if (file_stat->st_size > meta->largest_size) meta->largest_size = file_stat->st_size;

    if (has_dat_suffix(path) && file_stat->st_size > meta->largest_dat_size) {
        char *copy = strdup(path);
        if (!copy) {
            fprintf(stderr, "Out of memory while tracking the largest .dat file\n");
            exit(EXIT_FAILURE);
        }
        free(meta->largest_dat_path);
        meta->largest_dat_path = copy;
        meta->largest_dat_size = file_stat->st_size;
    }
    if (has_fits_suffix(path) && file_stat->st_size > meta->largest_fits_size) {
        char *copy = strdup(path);
        if (!copy) {
            fprintf(stderr, "Out of memory while tracking the largest .fits file\n");
            exit(EXIT_FAILURE);
        }
        free(meta->largest_fits_path);
        meta->largest_fits_path = copy;
        meta->largest_fits_size = file_stat->st_size;
    }
}

/* Size the NuFFT plan from at most one file per format: the largest .dat and the
 * largest .fits. Never inspect every file of a large batch - see AGENTS.md. */
static void finalize_batch_metadata(parameters *params, batch_meta_t *meta) {
    uint32_t max_newline_count = 0;
    double max_time_span = 0.0;

    if (meta->largest_dat_path) {
        uint32_t count = 0;
        double span = 0.0;
        inspect_dat_file(meta->largest_dat_path, &count, &span);
        if (count > max_newline_count) max_newline_count = count;
        if (span > max_time_span) max_time_span = span;
    }
    if (meta->largest_fits_path) {
        uint32_t count = 0;
        double span = 0.0;
        inspect_fits_file(meta->largest_fits_path, &count, &span);
        if (count > max_newline_count) max_newline_count = count;
        if (span > max_time_span) max_time_span = span;
    }

    if (max_newline_count == 0U) {
        fprintf(stderr, "No data rows found in the largest inspected .dat/.fits file of '%s'\n", params->target);
        free(meta->largest_dat_path);
        free(meta->largest_fits_path);
        meta->largest_dat_path = NULL;
        meta->largest_fits_path = NULL;
        exit(EXIT_FAILURE);
    }

    params->maxSize = (uint32_t)(meta->largest_size + 1);
    params->maxLen = max_newline_count;
    params->maxTimeSpan = max_time_span;

    size_t file_count = meta->file_count > 0U ? meta->file_count : 1U;
    params->avgLen = meta->total_bytes / ((uint64_t)DEFAULT_MEASUREMENT_SIZE * (uint64_t)file_count);
    params->maxFreqCount = estimate_frequency_count(params, params->maxTimeSpan);
    initialize_nufft_plan(params);

    free(meta->largest_dat_path);
    free(meta->largest_fits_path);
    meta->largest_dat_path = NULL;
    meta->largest_fits_path = NULL;
}

static void list_entry_skip_warning(size_t skipped, const char *entry, const char *reason) {
    if (skipped <= 10U) fprintf(stderr, "Warning: skipping list entry '%s' (%s)\n", entry, reason);
}

static void exit_with_list_cleanup(FILE *list, char *line, batch_meta_t *meta, const char *message) {
    fprintf(stderr, "%s\n", message);
    free(line);
    fclose(list);
    free(meta->largest_dat_path);
    free(meta->largest_fits_path);
    exit(EXIT_FAILURE);
}

/* Read a .tsv list: the first column is the photometry path, everything after the
 * first tab is ignored. Blank lines, '#' comments and the batch header line are
 * skipped so an ihsnpeaks .out.tsv can be fed back in directly. */
static void process_target_list(parameters *params) {
    FILE *list = fopen(params->target, "r");
    if (!list) {
        fprintf(stderr, "Failed to open input list '%s': %s\n", params->target, strerror(errno));
        exit(EXIT_FAILURE);
    }

    batch_meta_t meta;
    batch_meta_init(&meta);

    char *line = NULL;
    size_t line_capacity = 0;
    ssize_t line_length;
    size_t entry_count = 0U;
    size_t skipped_count = 0U;
    size_t absolute_count = 0U;
    bool first_line = true;

    while ((line_length = getline(&line, &line_capacity, list)) != -1) {
        while (line_length > 0 && (line[line_length - 1] == '\n' || line[line_length - 1] == '\r')) {
            line[--line_length] = '\0';
        }

        char *field = line;
        if (first_line) {
            if (line_length >= 3 && (unsigned char)field[0] == 0xEF && (unsigned char)field[1] == 0xBB && (unsigned char)field[2] == 0xBF) {
                field += 3;
            }
            first_line = false;
        }

        char *tab = strchr(field, '\t');
        if (tab) *tab = '\0';
        while (*field == ' ' || *field == '\t') field++;
        size_t field_len = strlen(field);
        while (field_len > 0U && (field[field_len - 1U] == ' ' || field[field_len - 1U] == '\t')) field[--field_len] = '\0';
        if (field_len == 0U || field[0] == '#' || strcmp(field, "Input_file") == 0) continue;

        entry_count++;

        bool bypassed = false;
        char *resolved = resolve_photometry_path(params, field, &bypassed);
        if (!resolved) exit_with_list_cleanup(list, line, &meta, "Out of memory while reading the input list");
        if (bypassed) absolute_count++;

        struct stat file_stat;
        if (stat(resolved, &file_stat) != 0) {
            int stat_errno = errno;
            char reason[160];
            snprintf(reason, sizeof(reason), "cannot stat: %s", strerror(stat_errno));
            skipped_count++;
            list_entry_skip_warning(skipped_count, field, reason);
            free(resolved);
            continue;
        }
        if (!S_ISREG(file_stat.st_mode)) {
            skipped_count++;
            list_entry_skip_warning(skipped_count, field, "not a regular file");
            free(resolved);
            continue;
        }
        if (!has_photometry_suffix(resolved)) {
            skipped_count++;
            list_entry_skip_warning(skipped_count, field, "unsupported format (expected .dat, .fit or .fits)");
            free(resolved);
            continue;
        }

        target_t target;
        target.label = strdup(field);
        if (!target.label) {
            free(resolved);
            exit_with_list_cleanup(list, line, &meta, "Out of memory while reading the input list");
        }
        target.path = resolved;
        batch_meta_add(&meta, resolved, &file_stat);
        kv_push(target_t, params->targets, target);
    }

    free(line);
    fclose(list);

    if (absolute_count > 0U && params->prefix && params->prefix[0] != '\0') {
        fprintf(stderr, "Warning: --prefix ignored for %zu absolute list entr%s\n", absolute_count, absolute_count == 1U ? "y" : "ies");
    }
    if (skipped_count > 0U) {
        fprintf(stderr, "Warning: skipped %zu of %zu list entr%s\n", skipped_count, entry_count, entry_count == 1U ? "y" : "ies");
    }
    if (meta.file_count == 0U) {
        free(meta.largest_dat_path);
        free(meta.largest_fits_path);
        fprintf(stderr, "No usable photometry files listed in '%s'\n", params->target);
        exit(EXIT_FAILURE);
    }

    finalize_batch_metadata(params, &meta);
}

/* Directory branch: labels use the target as typed; paths use the prefixed
 * directory, so --prefix only redirects where photometry is read from. */
static void process_target_directory(parameters *params, char *resolved_dir) {
    DIR *dir = opendir(resolved_dir);
    if (!dir) {
        fprintf(stderr, "Failed to open directory '%s': %s\n", resolved_dir, strerror(errno));
        free(resolved_dir);
        exit(EXIT_FAILURE);
    }

    batch_meta_t meta;
    batch_meta_init(&meta);

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        // d_type == DT_UNKNOWN is common on some filesystems; fall back to
        // stat(2) instead of silently skipping every candidate file.
        bool may_be_regular = (entry->d_type == DT_REG || entry->d_type == DT_UNKNOWN);
        if (!may_be_regular || !has_photometry_suffix(entry->d_name)) continue;

        char *label = join_directory_entry(params->target, entry->d_name);
        char *path = join_directory_entry(resolved_dir, entry->d_name);
        if (!label || !path) {
            free(label);
            free(path);
            closedir(dir);
            free(resolved_dir);
            fprintf(stderr, "Out of memory while scanning directory '%s'\n", params->target);
            exit(EXIT_FAILURE);
        }

        struct stat file_stat;
        if (stat(path, &file_stat) != 0) {
            int stat_errno = errno;
            fprintf(stderr, "Failed to stat '%s': %s\n", path, strerror(stat_errno));
            free(label);
            free(path);
            closedir(dir);
            free(resolved_dir);
            exit(EXIT_FAILURE);
        }
        if (!S_ISREG(file_stat.st_mode)) {
            free(label);
            free(path);
            continue;
        }

        target_t target = {.path = path, .label = label};
        batch_meta_add(&meta, path, &file_stat);
        kv_push(target_t, params->targets, target);
    }
    closedir(dir);
    free(resolved_dir);

    if (meta.file_count == 0U) {
        free(meta.largest_dat_path);
        free(meta.largest_fits_path);
        fprintf(stderr, "No regular .dat or .fits files found in '%s'\n", params->target);
        exit(EXIT_FAILURE);
    }

    finalize_batch_metadata(params, &meta);
}

static bool process_target_single(parameters *params, char *resolved) {
    target_t target;
    target.path = resolved;
    target.label = strdup(params->target);
    if (!target.label) {
        free(resolved);
        fprintf(stderr, "Out of memory\n");
        exit(EXIT_FAILURE);
    }

    struct stat file_stat;
    if (stat(resolved, &file_stat) != 0) {
        int stat_errno = errno;
        fprintf(stderr, "Failed to stat target '%s': %s\n", resolved, strerror(stat_errno));
        free(resolved);
        exit(EXIT_FAILURE);
    }
    params->maxSize = (uint32_t)(file_stat.st_size + 1);

    uint32_t line_count = 0U;
    double time_span = 0.0;
    if (has_dat_suffix(resolved)) {
        inspect_dat_file(resolved, &line_count, &time_span);
    } else if (has_fits_suffix(resolved)) {
        inspect_fits_file(resolved, &line_count, &time_span);
    } else {
        fprintf(stderr, "Unsupported input file format: '%s' (expected .dat, .fit, .fits or a .tsv list)\n", params->target);
        free(resolved);
        exit(EXIT_FAILURE);
    }
    if (line_count == 0U) {
        fprintf(stderr, "No data rows found in '%s'\n", params->target);
        free(resolved);
        exit(EXIT_FAILURE);
    }

    params->maxLen = line_count;
    params->maxTimeSpan = time_span;
    params->maxFreqCount = estimate_frequency_count(params, params->maxTimeSpan);
    initialize_nufft_plan(params);

    kv_push(target_t, params->targets, target);
    return true;
}

// Function to process the path and populate the target vector
static bool process_path(parameters *params) {
    kv_init(params->targets);
    params->maxLen = 0;
    params->maxTimeSpan = 0.0;
    params->avgLen = 0;
    params->isList = false;

    struct stat target_stat;
    bool target_exists = (stat(params->target, &target_stat) == 0);
    if (target_exists && S_ISREG(target_stat.st_mode) && has_tsv_suffix(params->target)) {
        params->isList = true;
        process_target_list(params);
        return false;
    }

    bool bypassed = false;
    char *resolved = resolve_photometry_path(params, params->target, &bypassed);
    if (!resolved) {
        fprintf(stderr, "Out of memory while resolving target '%s'\n", params->target);
        exit(EXIT_FAILURE);
    }
    if (bypassed) fprintf(stderr, "Warning: --prefix ignored for absolute target '%s'\n", params->target);

    struct stat resolved_stat;
    if (stat(resolved, &resolved_stat) != 0) {
        int stat_errno = errno;
        fprintf(stderr, "Failed to stat target '%s': %s\n", resolved, strerror(stat_errno));
        if (params->prefix && params->prefix[0] != '\0') {
            fprintf(stderr, "  (target '%s' with prefix '%s')\n", params->target, params->prefix);
        }
        free(resolved);
        exit(EXIT_FAILURE);
    }

    if (S_ISDIR(resolved_stat.st_mode)) {
        process_target_directory(params, resolved);
        return false;
    }
    return process_target_single(params, resolved);
}

// Function to free all memory allocated in the kvec and destroy it
static void free_targets(kvec_target_t *targets) {
    for (size_t i = 0; i < kv_size(*targets); i++) {
        free(kv_A(*targets, i).path);
        free(kv_A(*targets, i).label);
    }
    kv_destroy(*targets);
}

#endif  // METADATA_H
