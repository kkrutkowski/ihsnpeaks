#define DEFAULT_MEASUREMENT_SIZE 24 //typical for OGLE and MACHO .dat files. Change if needed.

#define MI_OVERRIDE 1
#include "../include/mimalloc/mimalloc.h"

#include <stdio.h>

#include "params.h"
#include "metadata.h"


int main(int argc, char *argv[]) {
    parameters params = read_parameters(argc, argv);

    // Initialize a kvec for target structs
    kvec_target_t targets;
    process_path(&params, &targets);

    // Output results
    //for (size_t i = 0; i < kv_size(targets); i++) {printf("File: %s\n", kv_A(targets, i).path);}
    //printf("File: %s\n", kv_A(targets, kv_size(targets)-1).path);

    // Output parsed parameters for verification
    printf("Parsed parameters:\n");
    printf("  Target: %s\n", params.target ? *params.target : "(none)");
    printf("  Maximum frequency: %.2f\n", params.fmax);
    printf("  Minimum frequency: %.2f\n", params.fmin);
    printf("  Oversampling Factor: %.2f\n", params.oversamplingFactor);
    printf("  Detection threshold: %.2f\n", params.treshold);
    printf("  Expected systemic variation: %.1e\n", params.epsilon);
    printf("  Npeaks: %d\n", params.npeaks);
    printf("  Nterms: %d\n", params.nterms);
    printf("  Spectrum: %s\n", params.spectrum ? "true" : "false");
    printf("  Is file: %s\n", params.isFile ? "true" : "false");

    // Free the allocated memory
    free_targets(&targets);
    free_parameters(&params);

    return 0;
}
