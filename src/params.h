#ifndef PARAMS_H
#define PARAMS_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

#include <klib/ketopt.h>

#include "utils/common.h"
#include "utils/readout.h"
#include "metadata.h"


// Function to initialize parameters with default values
static parameters init_parameters(int argc, char *argv[]) {
    parameters params = {0};

    params.target = (char *) malloc((strlen(argv[1]) + 1) * sizeof(char)); // +1 for the null terminator
    strcpy(params.target, argv[1]); // Copy the string

    params.fmax = atof(argv[2]);
    params.threshold = 10.0 * M_LN10;
    params.oversamplingFactor = 5.0;
    params.epsilon = 0.001;
    params.npeaks = 10;
    params.nterms = 3;
    params.defaultGridRatio = 32;
    params.gridRatio = 32;
    params.mode = 2;
    params.spectrum = false;

    return params;
}

void print_parameters(parameters *params) {
    printf("Parsed parameters:\n");
    printf("\tTarget: %s\n", params->target);
    printf("\tMaximum frequency: %.2f\n", params->fmax);
    printf("\tMinimum frequency: %.2f\n", params->fmin);
    printf("\tOversampling Factor: %.2f\n", params->oversamplingFactor);
    printf("\tDetection threshold: %.2f\n", params->threshold * M_LOG10E);
    printf("\tExpected systemic variation: %.1e\n", params->epsilon);
    printf("\tNpeaks: %d\n", params->npeaks);
    printf("\tNterms: %d\n", params->nterms);
    printf("\tSpectrum: %s\n", params->spectrum ? "true" : "false");
    printf("\tIs file: %s\n", params->isFile ? "true" : "false");
    printf("\tIDLE: %s\n", params->idle ? "true" : "false");
    printf("\tLargest file's length: %i\n", params->maxLen);
    printf("\tRead buffer size: %i\n", params->maxSize);
    printf("\tFFT grid length: %i\n", params->gridLen);
    printf("\tNumber of worker threads available: %i\n", params->jobs);
    if(!params->isFile){printf("\n");}
}

void alloc_buffers(parameters *params){
    params->buffers = calloc(params->nbuffers, sizeof(buffer_t*));
    for (int i = 0; i < params->nbuffers; i++) {params->buffers[i] = calloc(1, sizeof(buffer_t));}
}

// Free allocated memory for parameters
void free_parameters(parameters *params) {
    fftwf_destroy_plan(params->plan);
    free(params->target); // Free the allocated string
    free_targets(&params->targets);
    for (int i = 0; i < params->nbuffers; i++) {free_buffer(params->buffers[i]); free(params->buffers[i]);}
    free(params->buffers);
    if(params->outFile){free(params->outFile);}
}

void print_help(char** argv) {
    printf("Usage: %s target fmax [options]\n", argv[1]);
    printf("\n");
    printf("Positional arguments:\n");
    printf("  target                   Path containing .dat file(s) (file or directory)\n");
    printf("  fmax                     Upper bound of the frequency grid \n");
    printf("\n");
    printf("Options:\n");
    printf("  -o, --oversampling        Set expected number of frequencies per main lobe (default: 5.0)\n");
    printf("  -t, --threshold            Set the peak detection threshold (default: 10.0)\n");
    printf("  -f, --fmin                Lower bound of the frequency grid (default: 0.0)\n");
    printf("  -p, --peaks               Set the maximum number of peaks (default: 10)\n");
    printf("  -n, --terms               Set the number of harmonics used for computation (default: 3)\n");
    printf("  -j, --jobs                Limit of the number of worker threads used for computation (default: 0)\n");
    printf("  -m, --mode                Fraction of frequencies reevaluated with F-test (0-5, default:2)\n");
    printf("\n");
    printf("  -s, --spectrum            Print generated spectra into .tsv files (default: false)\n");
    printf("  -i, --idle                Use idle-type compute threads (default: false)\n");
    printf("\n"};
    printf("      --debug               Print parameters before the computation (default: false)\n");
    printf("\n"};
    printf("  -h, --help                Display this help message and exit\n");
    printf("Example:\n");
    printf("  %s /path/to/target.dat 10.0 -o 10.0 -t15.0 --peaks=5 -d --spectrum \n", argv[0]);
}

// Function to parse command-line arguments into a parameters struct
static parameters read_parameters(int argc, char *argv[]) {
    parameters params = init_parameters(argc, &argv[0]);

    // Define long options
    static ko_longopt_t longopts[] = {
         //impement the "generate mode" separately, precomputing the FFTW plans and saving them to /opt/ihnspeaks
        {"peaks", ko_required_argument, 'p'},
        {"terms", ko_required_argument, 'n'},
        {"threshold", ko_required_argument, 't'},
        {"fmin", ko_required_argument, 'f'},
        {"oversampling", ko_required_argument, 'o'},
        {"epsilon", ko_required_argument, 'e'},
        {"mode", ko_required_argument, 'm'},
        {"jobs", ko_required_argument, 'j'},
        {"spectrum", ko_no_argument, 's'},
        {"debug", ko_no_argument, 'd'},
        {"corrected", ko_no_argument, 'c'}, //apply the logarithmic correction, not fully implemented
        {"idle", ko_no_argument, 'i'},
        {"help", ko_no_argument, 'h'},
        // Negative cases, corresponding to long arguments only;
        {"pack", ko_no_argument, '\xfe'},
        {"strip", ko_no_argument, '\xfd'},
        {NULL, 0, 0}
    };

    // Initialize ketopt, skipping the first two positional arguments
    ketopt_t opt = KETOPT_INIT;
    opt.ind = 2; // Start parsing options from argv[2]

    int c;
    while ((c = ketopt(&opt, argc, argv, 1, "o:p:n:t:f:e:j:m:sdich\xfe\xfd", longopts)) != -1) {
        //printf("argument: %c", c);
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
                params.threshold = atof(opt.arg) * M_LN10;
                break;
            case 'f':
                params.fmin = atof(opt.arg);
                break;
            case 'e':
                params.epsilon = atof(opt.arg);
                break;
            case 'j':
                params.jobs = atoi(opt.arg);
                break;
            case 'm':
                params.mode = atoi(opt.arg);
                break;
            case 's':
                params.spectrum = true;
                break;
            case 'd':
                params.debug = true;
                break;
            case 'i':
                params.idle = true;
                break;
            case 'c':
                params.corrected = true;
                break;
            case 'h':
                print_help(argv); exit(0);
                break;
            case '\xfe':
                print_help(argv); exit(0);
                break;
            case '\xfd':
                print_help(argv); exit(0);
                break;
            case '?':
                fprintf(stderr, "Unknown option: -%c\n", opt.opt ? opt.opt : '?');
                break;
            case ':':
                fprintf(stderr, "Missing argument for -%c\n", opt.opt ? opt.opt : '?');
                break;
        }
    }
    params.isFile = process_path(&params);
    return params;
}

// Invert correct_ihs_res using the false position (regula falsi) method.
double correct_threshold(const double threshold, const int n) {
    const double tol = 1e-5;
    const int max_iter = 10;
    double low = 0, high = (threshold > 1.0 ? threshold : 1.0), x;
    while (correct_ihs_res(high, n) - threshold < 0) high *= 2;
    for (int i = 0; i < max_iter; i++) {
        double f_low = correct_ihs_res(low, n) - threshold;
        double f_high = correct_ihs_res(high, n) - threshold;
        x = low - f_low * (high - low) / (f_high - f_low);
        double f_x = correct_ihs_res(x, n) - threshold;
        if (fabs(f_x) < tol) break;
        if (f_x * f_low < 0) high = x; else low = x;
    }
    return x;
}


#endif // PARAMS_H
