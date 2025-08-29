#define DEFAULT_MEASUREMENT_SIZE 24
#define FFTW_MEASURE_THRESHOLD 19

#if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 202311L  // If C23 is not available
    #define constexpr const
#endif
#if (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 201112L)
    #define aligned_alloc aligned_alloc_fallback
#endif  /* aligned_alloc */

#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

static inline size_t bitceil(size_t n) {if(n <= 1){return 1;}
return (size_t) (2 << (int)ceil(log2(n)));}

/* C11-like aligned_alloc implementation for systems that don't already provide it.
   This function:
     - Rounds up the requested alignment to a power of two if needed,
     - Ensures the alignment is at least sizeof(void*) (required by posix_memalign),
     - Rounds up size to be a multiple of the (adjusted) alignment (as required by the C11 spec),
     - Allocates memory using posix_memalign.
   The returned pointer can be freed with free(). */
void *aligned_alloc_fallback(size_t alignment, size_t size) {
    void *ptr = NULL;
    if (alignment == 0) return NULL;
    if (alignment < sizeof(void*)) alignment = sizeof(void*);
    if ((alignment & (alignment - 1)) != 0) alignment = bitceil(alignment);
    if (size % alignment != 0) size = ((size + alignment - 1) / alignment) * alignment;
    if (posix_memalign(&ptr, alignment, size) != 0) return NULL;
    return ptr;
}

#include "params.h"
#include "metadata.h"
#include "utils/readout.h"
#include "process.h"

#include <mimalloc/mimalloc.h>
#include <klib/kthread.h>

int main(int argc, char *argv[]) {
    if (argc > 1) {
        if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "h") == 0) {print_help(argv); return 0;}
        else if (strcmp(argv[1], "generate") == 0 || strcmp(argv[1], "-h") == 0) {generate_plans(argv); return 0;}
    }
    if (argc < 3){print_help(argv); return 1;}
    if (access("/opt/ihsnpeaks/plans", F_OK) == 0) {fftwf_import_wisdom_from_filename("/opt/ihsnpeaks/plans");}
    else {printf("Warning: precomputed plans not found.\n Consider running 'sudo %s generate' for optimal performance.\n", argv[0]);}

    parameters params = read_parameters(argc, argv);
    int nThreads = sysconf(_SC_NPROCESSORS_ONLN);

    if (params.debug){print_parameters(&params);}
    if (params.mode < 5 && params.nterms > 1){params.threshold = correct_threshold(params.threshold, params.nterms);}

    // Output results
    //for (size_t i = 0; i < kv_size(params.targets); i++) {printf("File: %s\n", kv_A(params.targets, i).path);}
    //printf("File: %s\n", kv_A(params.targets, kv_size(params.targets)-1).path);

    if(kv_size(params.targets) == 1){
        //printf("Single file mode\n");
        buffer_t buffer = {0}; //initalize the pointers to NULL to avoid segfaults
        alloc_buffer(&buffer, &params); //, params.avgLen
            process_target(kv_A(params.targets, 0).path, &buffer, &params, false);
            //print_buffer(&buffer);
        //process the data here
        free_buffer(&buffer);
        free_parameters(&params);
    return 0;} //end the program's execution here if only one target is provided'

    //to avoid asigning more, than 1 thread per target
    else {
        // Extract the parent directory of params.target[0]
        char targetPath[256];
        strcpy(targetPath, params.target);  // Copy the target path to a modifiable buffer
        char *parentDir = dirname(targetPath); // Get the parent directory

        // Construct the full path for output.tsv
        char outputFilePath[256];
        snprintf(outputFilePath, sizeof(outputFilePath), "%s/output.tsv", parentDir);

        // Create the output.tsv file
        FILE *outputFile = fopen(outputFilePath, "w");
        if (outputFile == NULL) {perror("Failed to create the output file");return 1;}
        if (params.mode < 1){fprintf(outputFile, "Input_file\tbest_fit_frequency\t-log_10(FAP)\t\t[freq, -log_10(FAP)]\n");}
        else {fprintf(outputFile, "Input_file\tbest_fit_frequency\t-log_10(FAP)\t\t[freq, amp, R]\n");}
        fclose(outputFile);

        // Set the params.outFile pointer to the full path
        params.outFile = strdup(outputFilePath);
        printf("Output file's path: %s\n", outputFilePath); fflush(stdout);

        if (kv_size(params.targets) < nThreads) {nThreads = kv_size(params.targets);}
        if (params.jobs < nThreads && params.jobs > 0) {nThreads = params.jobs;}
        params.nbuffers = nThreads; alloc_buffers(&params);
        kt_forpool_t *pool = kt_forpool_init(nThreads, params.idle);

        //distribute processing of all of the targets on multiple threads
        kt_forpool(pool, process_targets, &params, kv_size(params.targets));
        for (int i = 0; i < params.nbuffers; i++) {fprint_buffer(params.buffers[i], &params);} //free the buffers before freeing

        // Kill the worker threads before exiting
        kt_forpool_destroy(pool);
        // Free the allocated memory
        free_parameters(&params);

    return 0;}
}
