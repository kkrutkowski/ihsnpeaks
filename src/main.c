#define DEFAULT_MEASUREMENT_SIZE 24
#include "../include/mimalloc/mimalloc.h"
#include <stdio.h>
#include <unistd.h>

#include "params.h"
#include "metadata.h"
#include "readout.h"
#include "../include/klib/kthread_batch.h"

int main(int argc, char *argv[]) {
    if (argc < 3){return 1;}
    parameters params = read_parameters(argc, argv);
    int nThreads = sysconf(_SC_NPROCESSORS_ONLN); void *thread_pool;

    // Initialize a kvec for target structs
    kvec_target_t targets; uint32_t maxLen; uint32_t maxSize;
    process_path(&params, &targets, &maxLen, &maxSize);
    print_parameters(&params);
    printf("  Largest file's length: %i\n", maxLen);
    printf("  Read buffer size: %i\n", maxSize);

    // Output results
    //for (size_t i = 0; i < kv_size(targets); i++) {printf("File: %s\n", kv_A(targets, i).path);}
    //printf("File: %s\n", kv_A(targets, kv_size(targets)-1).path);

    if(kv_size(targets) == 1){
        printf("Single file mode\n");
        buffer_t buffer = {0}; //initalize the pointers to NULL to avoid segfaults
        alloc_buffer(&buffer, maxLen, maxSize);
        read_dat(kv_A(targets, 0).path, &buffer);
        //process the data here
        free_buffer(&buffer); //segfaults
        free_targets(&targets);
        free_parameters(&params);
    return 0;} //end the program's execution here if only one target is provided'



    //to avoid asigning more, than 1 thread per target
    else {if (kv_size(targets) < nThreads) {nThreads = kv_size(targets);}}
    kt_forpool_t *fp = kt_forpool_init(nThreads);




    // Output parsed parameters for verification

    // Kill the worker threads before exiting
    kt_forpool_destroy(fp);
    // Free the allocated memory
    free_targets(&targets);
    free_parameters(&params);

    return 0;
}
