#ifndef PARAMS_H
#define PARAMS_H

//#include "../include/fast_convert.h"
#include "../include/klib/ketopt.h"

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
    int bufferSize;
    bool isFile;
    bool spectrum;
    bool debug;
} parameters;

// Function to initialize parameters with default values
static parameters init_parameters(int argc, char *argv[]) {
    parameters params;

    params.target = (char **) malloc(sizeof(char*)); // Allocate memory for the target string and copy the value from argv[1]
    params.target[0] = (char *) malloc((strlen(argv[1]) + 1) * sizeof(char)); // +1 for the null terminator
    strcpy(params.target[0], argv[1]); // Copy the string

    params.fmax = atof(argv[2]);
    params.fmin = 0.0;
    params.treshold = 4.0;
    params.oversamplingFactor = 3.0;
    params.epsilon = 0.001;
    params.npeaks = 10;
    params.nterms = 3;
    params.bufferSize = -1;
    params.spectrum = false;
    params.debug = false;

    return params;
}

// Free allocated memory for parameters
void free_parameters(parameters *params) {
    free(params->target[0]); // Free the allocated string
    free(params->target);     // Free the array of pointers
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
        {"bsize", ko_required_argument, 'b'},
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
    while ((c = ketopt(&opt, argc, argv, 1, "o:p:n:t:b:f:e:sd", longopts)) >= 0) {
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
            case 'b':
                params.bufferSize = atoi(opt.arg);
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

    return params;
}

#endif // PARAMS_H
