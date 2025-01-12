#define DEFAULT_MEASUREMENT_SIZE 24
#include "../include/mimalloc/mimalloc.h"
#include <stdio.h>
#include <unistd.h>
#include <libgen.h>

#include "params.h"
#include "metadata.h"
#include "utils/readout.h"
#include "process.h"

#include "../include/klib/kthread.h"

int main(int argc, char *argv[]) {
    if (argc < 3){return 1;}
    parameters params = read_parameters(argc, argv);
    int nThreads = sysconf(_SC_NPROCESSORS_ONLN); void *thread_pool;

    if (params.debug){print_parameters(&params);}

    // Output results
    //for (size_t i = 0; i < kv_size(params.targets); i++) {printf("File: %s\n", kv_A(params.targets, i).path);}
    //printf("File: %s\n", kv_A(params.targets, kv_size(params.targets)-1).path);

    if(kv_size(params.targets) == 1){
        //printf("Single file mode\n");
        buffer_t buffer = {0}; //initalize the pointers to NULL to avoid segfaults
        alloc_buffer(&buffer, params.nterms, params.maxLen, params.maxSize, params.gridLen, params.npeaks); //, params.avgLen
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
        fclose(outputFile);

        // Set the params.outFile pointer to the full path
        params.outFile = strdup(outputFilePath);

        if (kv_size(params.targets) < nThreads) {nThreads = kv_size(params.targets);}
        params.nbuffers = nThreads; alloc_buffers(&params);
        kt_forpool_t *fp = kt_forpool_init(nThreads);


        // Output parsed parameters for verification

        // Kill the worker threads before exiting
        kt_forpool_destroy(fp);
        // Free the allocated memory
        free_parameters(&params);

    return 0;}
}
