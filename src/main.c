#define DEFAULT_MEASUREMENT_SIZE 24
#define IHSNPEAKS_VERSION "v1.1.0-preview"

#include <klib/kthread.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include "utils/compat.h"

#ifndef HAS_MIMALLOC
#    define HAS_MIMALLOC 0
#endif

#if HAS_MIMALLOC
#    include <mimalloc/mimalloc.h>
#endif

#include "metadata.h"
#include "params.h"
#include "process.h"
#include "utils/readout.h"

#include "hwloc_topo.h"
#include "pool.h"

static void build_batch_output_path(const char *target, char *out_file, size_t out_file_size) {
    char target_path[PATH_MAX];
    strncpy(target_path, target, sizeof(target_path) - 1U);
    target_path[sizeof(target_path) - 1U] = '\0';

    size_t len = strlen(target_path);
    while (len > 1U && target_path[len - 1U] == '/') {
        target_path[--len] = '\0';
    }

    char *slash = strrchr(target_path, '/');
    char *dot = strrchr(target_path, '.');
    if (dot && (!slash || dot > slash)) {
        *dot = '\0';
    }

    snprintf(out_file, out_file_size, "%s.tsv", target_path);
}

int main(int argc, char *argv[]) {
    if (argc > 1) {
        if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "h") == 0) {
            print_help(argv);
            return 0;
        } else if (strcmp(argv[1], "generate") == 0 || strcmp(argv[1], "--generate") == 0) {
            // Early-exit topology generation (no target/fmax required)
            int nThreads = sysconf(_SC_NPROCESSORS_ONLN);
            if (nThreads < 1) nThreads = 1;
            topo_config_t *topo = probe_topology(nThreads);
            if (topo) {
                save_topology_config(topo);  // prints summary internally
                free_topology(topo);
            } else {
                fprintf(stderr, "Failed to probe hardware topology\n");
                return 1;
            }
            return 0;
        }
    }
    if (argc < 3) {
        print_help(argv);
        return 1;
    }
    parameters params = read_parameters(argc, argv);
    int nThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nThreads < 1) nThreads = 1;
    if (params.jobs < nThreads && params.jobs > 0) {
        nThreads = params.jobs;
    }

    // Probe hardware topology (L3 domains, cores, SMT).
    // Auto-saves ~/.ihsnpeaks/hwloc_config on first run.
    topo_config_t *topo = load_or_probe_topology(params.generate, nThreads);
    if (!topo) {
        fprintf(stderr, "Failed to probe hardware topology via hwloc\n");
        return 1;
    }
    if (params.generate) {
        free_topology(topo);
        free_parameters(&params);
        return 0;
    }
    nThreads = topo->total_workers;
    if (nThreads < 1) nThreads = 1;

    if (params.debug) {
        print_topology_summary(topo);
        print_parameters(&params);
    }
    // Output results
    // for (size_t i = 0; i < kv_size(params.targets); i++) {printf("File: %s\n", kv_A(params.targets, i).path);}
    // printf("File: %s\n", kv_A(params.targets, kv_size(params.targets)-1).path);

    if (kv_size(params.targets) == 1) {
        // printf("Single file mode\n");
        buffer_t buffer = {0};           // initalize the pointers to NULL to avoid segfaults
        alloc_buffer(&buffer, &params);  //, params.avgLen
        kt_forpool_t *directPool = NULL;
        if (mode_uses_direct_eval_grid(params.mode) && nThreads > 1) {
            directPool = kt_forpool_init(nThreads, params.idle);
        }
        process_target(kv_A(params.targets, 0).path, &buffer, &params, false, directPool);
        // print_buffer(&buffer);
        // process the data here
        profile_report(&params);
        if (directPool) kt_forpool_destroy(directPool);
        free_buffer(&buffer);
        free_topology(topo);
        free_parameters(&params);
        return 0;
    }  // end the program's execution here if only one target is provided'

    // to avoid asigning more, than 1 thread per target
    else {
        char outputFilePath[PATH_MAX];
        build_batch_output_path(params.target, outputFilePath, sizeof(outputFilePath));

        // Create the output file
        FILE *outputFile = fopen(outputFilePath, "w");
        if (outputFile == NULL) {
            perror("Failed to create the output file");
            return 1;
        }
        const char *coordinateLabel = params.outputPeriod ? "period" : "freq";
        const char *bestFitLabel = params.outputPeriod ? "best_fit_period" : "best_fit_frequency";
        const eval_method_t *evalMethod = eval_method_for_params(&params);
        const char *powerLabel = evalMethod->peak_power_label;
        if (periodogram_uses_aov(params.periodogramMethod) && params.mode < 1) {
            fprintf(outputFile, "Input_file\t%s\t%s\t\t[%s, %s, R2]\n", bestFitLabel, powerLabel, coordinateLabel, powerLabel);
        } else if (params.mode < 1) {
            fprintf(outputFile, "Input_file\t%s\t%s\t\t[%s, %s]\n", bestFitLabel, powerLabel, coordinateLabel, powerLabel);
        } else {
            if (evalMethod->width_label) {
                fprintf(outputFile, "Input_file\t%s\t%s\t\t[%s, amp, %s, %s]\n", bestFitLabel, powerLabel, coordinateLabel, evalMethod->stat_label,
                        evalMethod->width_label);
            } else {
                fprintf(outputFile, "Input_file\t%s\t%s\t\t[%s, amp, %s]\n", bestFitLabel, powerLabel, coordinateLabel, evalMethod->stat_label);
            }
        }
        fclose(outputFile);

        // Set the params.outFile pointer to the full path
        params.outFile = strdup(outputFilePath);
        printf("Output file's path: %s\n", outputFilePath);
        fflush(stdout);

        size_t targetCount = kv_size(params.targets);
        bool splitDirectBatch = false;
        if (mode_uses_direct_eval_grid(params.mode) && nThreads > 0) {
            size_t filesPerWorker = (targetCount + (size_t)nThreads - 1U) / (size_t)nThreads;
            splitDirectBatch = filesPerWorker < 5U;
        }

        if (splitDirectBatch) {
            params.nbuffers = 1;
            alloc_buffers(&params);
            kt_forpool_t *directPool = NULL;
            if (nThreads > 1) {
                directPool = kt_forpool_init(nThreads, params.idle);
            }
            int permile = (int)targetCount / 1000;
            if (permile == 0) permile = 1;
            if (!params.buffers[0]->allocated) {
                alloc_buffer(params.buffers[0], &params);
            }
            for (size_t i = 0; i < targetCount; ++i) {
                int current = (int)i + 1;
                if (current % permile == 0 || (size_t)current == targetCount) {
                    float progress = (float)current * 100.0f / (float)targetCount;
                    printf("Computation in progress: %.1f%% complete\r", progress);
                    fflush(stdout);
                }
                process_target(kv_A(params.targets, i).path, params.buffers[0], &params, true, directPool);
                fprint_buffer(params.buffers[0], &params);
            }
            profile_report(&params);
            if (directPool) kt_forpool_destroy(directPool);
            free_topology(topo);
            free_parameters(&params);
            return 0;
        }

        if (targetCount < (size_t)topo->total_workers) {
            // Fewer files than total workers: use direct single-pool approach
            // to avoid empty pools. Re-probe with capped thread count.
            free_topology(topo);
            topo = probe_topology((int)(targetCount > INT32_MAX ? INT32_MAX : (long)targetCount));
            if (!topo) {
                fprintf(stderr, "Failed to re-probe topology for small batch\n");
                return 1;
            }
            nThreads = topo->total_workers;
            if (nThreads < 1) nThreads = 1;
        }

        // L3-aware dispatch: one pool per L3 domain, chunked file distribution
        l3_dispatcher_t *disp = l3_dispatcher_create(&params, topo);

        // Distribute processing of all targets across L3-local pools
        l3_dispatcher_run(disp);

        for (int i = 0; i < params.nbuffers; i++) {
            fprint_buffer(params.buffers[i], &params);
        }

        profile_report(&params);

        l3_dispatcher_destroy(disp);
        free_topology(topo);
        free_parameters(&params);

        return 0;
    }
}
