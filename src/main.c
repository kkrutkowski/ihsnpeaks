#define MI_OVERRIDE 1

#include "params.h"
#include "readout.h"

#include <stdio.h>
#include "../include/mimalloc/mimalloc.h"

int main(int argc, char *argv[]) {
    parameters params = read_parameters(argc, argv);


    // Output parsed parameters for verification
    printf("Parsed parameters:\n");
    printf("  Target: %s\n", params.target ? *params.target : "(none)");
    printf("  Maximum frequency: %.2f\n", params.fmax);
    printf("  Minimum frequency: %.2f\n", params.fmin);
    printf("  Oversampling Factor: %.2f\n", params.oversamplingFactor);
    printf("  Detection treshold: %.2f\n", params.treshold);
    printf("  Spectrum: %s\n", params.spectrum ? "true" : "false");
    printf("  Npeaks: %d\n", params.npeaks);
    printf("  Nterms: %d\n", params.nterms);

    free_parameters(&params);

    return 0;
}
