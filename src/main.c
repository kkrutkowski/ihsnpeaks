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

    print_parameters(&params);

    // Output results
    //for (size_t i = 0; i < kv_size(params.targets); i++) {printf("File: %s\n", kv_A(params.targets, i).path);}
    //printf("File: %s\n", kv_A(params.targets, kv_size(params.targets)-1).path);

    if(kv_size(params.targets) == 1){
        printf("Single file mode\n");
        buffer_t buffer = {0}; //initalize the pointers to NULL to avoid segfaults
        //alloc_buffer(&buffer, params.maxLen, params.maxSize);
        //read_dat(kv_A(params.targets, 0).path, &buffer);
        //process the data here
        //free_buffer(&buffer);
        free_parameters(&params);
    return 0;} //end the program's execution here if only one target is provided'



    //to avoid asigning more, than 1 thread per target
    else {if (kv_size(params.targets) < nThreads) {nThreads = kv_size(params.targets);}}
    kt_forpool_t *fp = kt_forpool_init(nThreads);


    // Output parsed parameters for verification

    // Kill the worker threads before exiting
    kt_forpool_destroy(fp);
    // Free the allocated memory
    free_parameters(&params);

    return 0;
}
