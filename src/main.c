#define DEFAULT_MEASUREMENT_SIZE 24
#define IHSNPEAKS_VERSION "v1.0.0-preview"

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

int main(int argc, char *argv[]) {
    if (argc > 1) {
        if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "h") == 0) {
            print_help(argv);
            return 0;
        } else if (strcmp(argv[1], "generate") == 0 || strcmp(argv[1], "-h") == 0) {
            generate_plans(argv);
            return 0;
        }
    }
    if (argc < 3) {
        print_help(argv);
        return 1;
    }
    parameters params = read_parameters(argc, argv);
    int nThreads = sysconf(_SC_NPROCESSORS_ONLN);

    if (params.debug) {
        print_parameters(&params);
    }
    // Output results
    // for (size_t i = 0; i < kv_size(params.targets); i++) {printf("File: %s\n", kv_A(params.targets, i).path);}
    // printf("File: %s\n", kv_A(params.targets, kv_size(params.targets)-1).path);

    if (kv_size(params.targets) == 1) {
        // printf("Single file mode\n");
        buffer_t buffer = {0};           // initalize the pointers to NULL to avoid segfaults
        alloc_buffer(&buffer, &params);  //, params.avgLen
        process_target(kv_A(params.targets, 0).path, &buffer, &params, false);
        // print_buffer(&buffer);
        // process the data here
        profile_report(&params);
        free_buffer(&buffer);
        free_parameters(&params);
        return 0;
    }  // end the program's execution here if only one target is provided'

    // to avoid asigning more, than 1 thread per target
    else {
        // Extract the parent directory of params.target[0]
        char targetPath[256];
        strcpy(targetPath, params.target);      // Copy the target path to a modifiable buffer
        char *parentDir = dirname(targetPath);  // Get the parent directory

        // Construct the full path for output.tsv
        char outputFilePath[256];
        snprintf(outputFilePath, sizeof(outputFilePath), "%s/output.tsv", parentDir);

        // Create the output.tsv file
        FILE *outputFile = fopen(outputFilePath, "w");
        if (outputFile == NULL) {
            perror("Failed to create the output file");
            return 1;
        }
        const char *coordinateLabel = params.outputPeriod ? "period" : "freq";
        const char *bestFitLabel = params.outputPeriod ? "best_fit_period" : "best_fit_frequency";
        if (periodogram_uses_aov(params.periodogramMethod) && params.mode < 1) {
            fprintf(outputFile, "Input_file\t%s\t-log_10(FAP)\t\t[%s, -log_10(FAP), R2]\n", bestFitLabel, coordinateLabel);
        } else if (params.mode < 1) {
            fprintf(outputFile, "Input_file\t%s\t-log_10(FAP)\t\t[%s, -log_10(FAP)]\n", bestFitLabel, coordinateLabel);
        } else {
            fprintf(outputFile, "Input_file\t%s\t-log_10(FAP)\t\t[%s, amp, %s]\n", bestFitLabel, coordinateLabel, gb_stat_label(params.gbEvalMode));
        }
        fclose(outputFile);

        // Set the params.outFile pointer to the full path
        params.outFile = strdup(outputFilePath);
        printf("Output file's path: %s\n", outputFilePath);
        fflush(stdout);

        if (kv_size(params.targets) < nThreads) {
            nThreads = kv_size(params.targets);
        }
        if (params.jobs < nThreads && params.jobs > 0) {
            nThreads = params.jobs;
        }
        params.nbuffers = nThreads;
        alloc_buffers(&params);
        kt_forpool_t *pool = kt_forpool_init(nThreads, params.idle);

        // distribute processing of all of the targets on multiple threads
        kt_forpool(pool, process_targets, &params, kv_size(params.targets));
        for (int i = 0; i < params.nbuffers; i++) {
            fprint_buffer(params.buffers[i], &params);
        }  // free the buffers before freeing

        profile_report(&params);

        // Kill the worker threads before exiting
        kt_forpool_destroy(pool);
        // Free the allocated memory
        free_parameters(&params);

        return 0;
    }
}
