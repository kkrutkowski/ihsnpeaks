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
    kvec_target_t targets; uint32_t maxLen;
    process_path(&params, &targets, &maxLen);
    print_parameters(&params);
    printf("  Largest file's size: %i\n", maxLen);

    // Output results
    //for (size_t i = 0; i < kv_size(targets); i++) {printf("File: %s\n", kv_A(targets, i).path);}
    //printf("File: %s\n", kv_A(targets, kv_size(targets)-1).path);

    if(kv_size(targets) == 1){
        printf("Single file mode\n");
        //process the data here
        free_targets(&targets);
        free_parameters(&params);
    return EXIT_SUCCESS;} //end the program's execution here if only one target is provided'



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
