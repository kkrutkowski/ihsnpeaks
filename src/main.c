#define DEFAULT_MEASUREMENT_SIZE 24 //typical for OGLE and MACHO .dat files. Change if needed.

#define MI_OVERRIDE 1
#include "../include/mimalloc/mimalloc.h"

#include <stdio.h>

#include "params.h"
#include "metadata.h"


int main(int argc, char *argv[]) {
    if (argc < 3){return 1;}
    parameters params = read_parameters(argc, argv);

    print_parameters(&params);

    // Initialize a kvec for target structs
    kvec_target_t targets;
    process_path(&params, &targets);

    // Output results
    //for (size_t i = 0; i < kv_size(targets); i++) {printf("File: %s\n", kv_A(targets, i).path);}
    //printf("File: %s\n", kv_A(targets, kv_size(targets)-1).path);

    // Output parsed parameters for verification


    // Free the allocated memory
    free_targets(&targets);
    free_parameters(&params);

    return 0;
}
