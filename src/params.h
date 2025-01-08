#ifndef PARAMS_H
#define PARAMS_H

//#include "../include/fast_convert.h"
#include "../include/klib/ketopt.h"

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "metadata.h"
#include "readout.h"

typedef struct {
    char**  target;
    float** buffer;
    float fmin;
    float fmax;
    float treshold;
    float oversamplingFactor;
    float epsilon;
    int npeaks;
    int nterms;
    bool isFile;
    bool spectrum;
    bool debug;

    fftwf_plan plan;

    //Variables used to estimate optimal hyperparameters from metadata
    //uint32_t blockSize;       //size of the FFT plan applied
    //uint32_t bufferSize;      //size of the FFT buffer
    uint32_t maxSize;         //number of bytes in the longest time series of the processed batch
    uint32_t maxLen;          //number of measurements in the longest time series of the processed batch
    uint32_t gridLen;         //length of the FFT grid used in transform
    uint64_t avgLen;          //average number of measurements per time series in the batch

    kvec_target_t targets;

    buffer_t** buffers;
} parameters;

// Function to initialize parameters with default values
static parameters init_parameters(int argc, char *argv[]) {
    parameters params = {0};

    params.target = (char **) malloc(sizeof(char*)); // Allocate memory for the target string and copy the value from argv[1]
    params.target[0] = (char *) malloc((strlen(argv[1]) + 1) * sizeof(char)); // +1 for the null terminator
    strcpy(params.target[0], argv[1]); // Copy the string

    params.fmax = atof(argv[2]);
    params.fmin = 0.0;
    params.treshold = 4.0;
    params.oversamplingFactor = 5.0;
    params.epsilon = 0.001;
    params.npeaks = 10;
    params.nterms = 3;
    params.spectrum = false;
    params.debug = false;

    return params;
}

void print_parameters(parameters *params) {
    printf("Parsed parameters:\n");
    printf("  Target: %s\n", params->target ? *params->target : "(none)");
    printf("  Maximum frequency: %.2f\n", params->fmax);
    printf("  Minimum frequency: %.2f\n", params->fmin);
    printf("  Oversampling Factor: %.2f\n", params->oversamplingFactor);
    printf("  Detection threshold: %.2f\n", params->treshold);
    printf("  Expected systemic variation: %.1e\n", params->epsilon);
    printf("  Npeaks: %d\n", params->npeaks);
    printf("  Nterms: %d\n", params->nterms);
    printf("  Spectrum: %s\n", params->spectrum ? "true" : "false");
    printf("  Is file: %s\n", params->isFile ? "true" : "false");
    printf("  Largest file's length: %i\n", params->maxLen);
    printf("  Read buffer size: %i\n", params->maxSize);
}

// Free allocated memory for parameters
void free_parameters(parameters *params) {
    fftwf_destroy_plan(params->plan);
    free(params->target[0]); // Free the allocated string
    free(params->target);     // Free the array of pointers
    free_targets(&params->targets);
}


// Function to parse command-line arguments into a parameters struct
static parameters read_parameters(int argc, char *argv[]) {
    parameters params = init_parameters(argc, &argv[0]);
    //params.target = &argv[1];


    // Define long options
    static ko_longopt_t longopts[] = {
        {"peaks", ko_required_argument, 'p'},
        {"terms", ko_required_argument, 'n'},
        {"treshold", ko_required_argument, 't'},
        {"msize", ko_required_argument, 'm'},
        {"fmin", ko_required_argument, 'f'},
        {"oversampling", ko_required_argument, 'o'},
        {"epsilon", ko_required_argument, 'e'},
        {"spectrum", ko_no_argument, 's'},
        {"debug", ko_no_argument, 'd'},
        {NULL, 0, 0}
    };

    // Initialize ketopt, skipping the first two positional arguments
    ketopt_t opt = KETOPT_INIT;
    opt.ind = 3; // Start parsing options from argv[3]

    int c;
    while ((c = ketopt(&opt, argc, argv, 1, "o:p:n:t:m:f:e:sd", longopts)) >= 0) {
        switch (c) {
            case 'o':
                params.oversamplingFactor = atof(opt.arg);
                break;
            case 'p':
                params.npeaks = atoi(opt.arg);
                break;
            case 'n':
                params.nterms = atoi(opt.arg);
                break;
            case 't':
                params.treshold = atof(opt.arg);
                break;
            case 'f':
                params.fmin = atof(opt.arg);
                break;
            case 'e':
                params.epsilon = atof(opt.arg);
                break;
            case 's':
                params.spectrum = true;
                break;
            case 'd':
                params.debug = true;
                break;
            case '?':
                fprintf(stderr, "Unknown option: -%c\n", opt.opt ? opt.opt : '?');
                break;
            case ':':
                fprintf(stderr, "Missing argument for -%c\n", opt.opt ? opt.opt : '?');
                break;
        }
    }

    params.isFile = process_path(params.target[0], &params.targets, &params.maxLen, &params.maxSize, &params.avgLen, &params.gridLen, &params.plan);

    return params;}

#endif // PARAMS_H
