#ifndef PARAMS_H
#define PARAMS_H

#include <errno.h>
#include <klib/ketopt.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "metadata.h"
#include "utils/common.h"
#include "utils/readout.h"

// Function to initialize parameters with default values
static parameters init_parameters(int argc, char *argv[]) {
    parameters params = {0};

    params.target = (char *)malloc((strlen(argv[1]) + 1) * sizeof(char));  // +1 for the null terminator
    strcpy(params.target, argv[1]);                                        // Copy the string

    params.fmax = atof(argv[2]);
    params.threshold = 10.0 * M_LN10;
    params.oversamplingFactor = 5.0;
    params.epsilon = 0.001;
    params.gbAlpha = 0.025f;
    params.npeaks = 10;
    params.nterms = 3;
    params.defaultGridRatio = 128;
    params.gridRatio = 128;
    params.gridMode = NUFFT1_PSWF43;
    params.mode = 2;
    params.r2_threshold = 0.05;
    params.periodogramMethod = PERIODOGRAM_IHS;
    params.gbEvalMode = GB_EVAL_GBLS;

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
    printf("\tPeriodogram method: %s\n", periodogram_method_name(params->periodogramMethod));
    printf("\tGaussian blur evaluation: %s\n", gb_eval_name(params->gbEvalMode));
    printf("\tGaussian blur alpha: %.6g\n", params->gbAlpha);
    printf("\tNpeaks: %d\n", params->npeaks);
    printf("\tNterms: %d\n", params->nterms);
    printf("\tSpectrum: %s\n", params->spectrum ? "true" : "false");
    printf("\tPrewhiten: %s\n", params->prewhiten ? "true" : "false");
    printf("\tIs file: %s\n", params->isFile ? "true" : "false");
    printf("\tIDLE: %s\n", params->idle ? "true" : "false");
    printf("\tLargest file's length: %i\n", params->maxLen);
    printf("\tRead buffer size: %i\n", params->maxSize);
    printf("\tNuFFT grid: pswf%i\n", params->gridMode);
    printf("\tNuFFT plan length: %i\n", params->gridLen);
    printf("\tNuFFT output block length: %i\n", params->outputLen);
    printf("\tNuFFT twiddle ladder levels: %i\n", params->ladderLevels);
    printf("\tMaximum target frequencies: %i\n", params->maxFreqCount);
    printf("\tNumber of worker threads available: %i\n", params->jobs);
    if (!params->isFile) {
        printf("\n");
    }
}

void alloc_buffers(parameters *params) {
    params->buffers = calloc(params->nbuffers, sizeof(buffer_t *));
    for (int i = 0; i < params->nbuffers; i++) {
        params->buffers[i] = calloc(1, sizeof(buffer_t));
    }
}

// Free allocated memory for parameters
void free_parameters(parameters *params) {
    nufft1_free_plan(params->nufftPlan);
    free(params->nufftTwiddleReal);
    free(params->nufftTwiddleImag);
    free(params->target);  // Free the allocated string
    free_targets(&params->targets);
    for (int i = 0; i < params->nbuffers; i++) {
        free_buffer(params->buffers[i]);
        free(params->buffers[i]);
    }
    free(params->buffers);
    if (params->outFile) {
        free(params->outFile);
    }
}

void print_help(char **argv) {
    printf("Usage: %s target fmax [options]\n", argv[0]);
    printf("\n");
    printf("Positional arguments:\n");
    printf("  target                   Path containing .dat file(s) (file or directory)\n");
    printf("  fmax                     Upper bound of the frequency grid \n");
    printf("\n");
    printf("Options:\n");
    printf("  -o, --oversampling        Set expected number of frequencies per main lobe (default: 5.0)\n");
    printf("  -t, --threshold            Set the peak detection threshold (default: 10.0)\n");
    printf("  -f, --fmin                Lower bound of the frequency grid (default: 0.0)\n");
    printf("  -n, --peaks               Set the maximum number of peaks (default: 10)\n");
    printf("  -d, --terms               Set the number of harmonics used for computation (default: 3)\n");
    printf("  -j, --jobs                Limit of the number of worker threads used for computation (default: 0)\n");
    printf("  -m, --mode                Peak evaluation mode (0-6, default:2):\n");
    printf("                            0 NuFFT-based grid only\n");
    printf("                            1 reevaluate up to n strongest above-threshold grid peaks\n");
    printf("                            2 like 1, plus binary search maximum refinement\n");
    printf("                            3 evaluate all local peaks, then binary-search the n strongest\n");
    printf("                            4 evaluate and binary-search every local peak\n");
    printf("                            5 like 2, but use a dense direct evaluation grid\n");
    printf("                            6 like 4, but dense direct evaluation replaces the NuFFT grid\n");
    printf("  -g, --g, --grid           Periodogram method: ihs (default) or aov/aovmh/aobmhw/chi/chi2/fastchi2\n");
    printf("  -e, --eval, --evaluate    Gaussian blur evaluation: gbls|gbl or gbas|gba, optionally [alpha in 0..1] (default: gbls[0.025])\n");
    printf("      --epsilon             Set expected systemic variation (default: 0.001)\n");
    printf("      --nfft, --nufft, --nufft1\n");
    printf("                            NuFFT backend: 43|pswf43 or 21|pswf21 (default: pswf43)\n");
    printf("\n");
    printf("  -s, --spectrum            Print generated spectra into .tsv files (default: false)\n");
    printf("  -i, --idle                Use idle-type compute threads (default: false)\n");
    printf("\n");
    printf("      --debug               Print parameters before the computation (default: false)\n");
    printf("      --prewhiten           Attenuate detected variability modes (default: false)\n");
    printf("\n");
    printf("  -h, --help                Display this help message and exit\n");
    printf("Example:\n");
    printf("  %s /path/to/target.dat 10.0 -o 10.0 -t15.0 --peaks=5 --debug --s \n", argv[0]);
}

static bool parse_nufft_mode(const char *arg, nufft1_mode *mode) {
    if (strcmp(arg, "43") == 0 || strcmp(arg, "pswf43") == 0) {
        *mode = NUFFT1_PSWF43;
        return true;
    }
    if (strcmp(arg, "21") == 0 || strcmp(arg, "pswf21") == 0) {
        *mode = NUFFT1_PSWF21;
        return true;
    }
    return false;
}

static bool parse_periodogram_method(const char *arg, periodogram_method *method) {
    if (strcmp(arg, "ihs") == 0) {
        *method = PERIODOGRAM_IHS;
        return true;
    }
    if (strcmp(arg, "aov") == 0) {
        *method = PERIODOGRAM_AOV;
        return true;
    }
    if (strcmp(arg, "aovmh") == 0) {
        *method = PERIODOGRAM_AOVMH;
        return true;
    }
    if (strcmp(arg, "aobmhw") == 0) {
        *method = PERIODOGRAM_AOBMHW;
        return true;
    }
    if (strcmp(arg, "chi") == 0) {
        *method = PERIODOGRAM_CHI;
        return true;
    }
    if (strcmp(arg, "chi2") == 0) {
        *method = PERIODOGRAM_CHI2;
        return true;
    }
    if (strcmp(arg, "fastchi2") == 0) {
        *method = PERIODOGRAM_FASTCHI2;
        return true;
    }
    return false;
}

static bool parse_gb_eval(const char *arg, gb_eval_mode *mode, float *gbAlpha) {
    const char *open = strchr(arg, '[');
    const char *name_end = open ? open : arg + strlen(arg);
    size_t name_len = (size_t)(name_end - arg);
    char name[16];
    if (name_len == 0 || name_len >= sizeof(name)) return false;
    memcpy(name, arg, name_len);
    name[name_len] = '\0';

    if (open) {
        const char *close = strchr(open + 1, ']');
        if (!close || close[1] != '\0') return false;
        errno = 0;
        char *parse_end = NULL;
        float parsed_alpha = strtof(open + 1, &parse_end);
        if (errno != 0 || parse_end != close || !isfinite(parsed_alpha)) return false;
        if (parsed_alpha < 0.0f) {
            fprintf(stderr, "Warning: Gaussian blur alpha %.6g is below 0.0; clamping to 0.0.\n", parsed_alpha);
            parsed_alpha = 0.0f;
        }
        if (parsed_alpha > 1.0f) {
            fprintf(stderr, "Warning: Gaussian blur alpha %.6g is above 1.0; clamping to 1.0.\n", parsed_alpha);
            parsed_alpha = 1.0f;
        }
        *gbAlpha = parsed_alpha;
    }

    if (strcmp(name, "gbls") == 0 || strcmp(name, "gbl") == 0) {
        *mode = GB_EVAL_GBLS;
        return true;
    }
    if (strcmp(name, "gbas") == 0 || strcmp(name, "gba") == 0) {
        *mode = GB_EVAL_GBAS;
        return true;
    }
    return false;
}

// Function to parse command-line arguments into a parameters struct
static parameters read_parameters(int argc, char *argv[]) {
    parameters params = init_parameters(argc, &argv[0]);

    enum { OPT_EPSILON = 0xF9, OPT_NUFFT = 0xFA, OPT_PREWHITEN = 0xFB, OPT_DEBUG = 0xFC, OPT_STRIP = 0xFD, OPT_PACK = 0xFE };

    // Define long options
    static ko_longopt_t longopts[] = {// impement the "generate mode" separately, precomputing the FFTW plans and saving them to /opt/ihnspeaks
                                      {"peaks", ko_required_argument, 'n'},
                                      {"terms", ko_required_argument, 'd'},
                                      {"threshold", ko_required_argument, 't'},
                                      {"fmin", ko_required_argument, 'f'},
                                      {"oversampling", ko_required_argument, 'o'},
                                      {"epsilon", ko_required_argument, OPT_EPSILON},
                                      {"eval", ko_required_argument, 'e'},
                                      {"evaluate", ko_required_argument, 'e'},
                                      {"mode", ko_required_argument, 'm'},
                                      {"g", ko_required_argument, 'g'},
                                      {"grid", ko_required_argument, 'g'},
                                      {"nfft", ko_required_argument, OPT_NUFFT},
                                      {"nufft", ko_required_argument, OPT_NUFFT},
                                      {"nufft1", ko_required_argument, OPT_NUFFT},
                                      {"jobs", ko_required_argument, 'j'},
                                      {"spectrum", ko_no_argument, 's'},
                                      {"corrected", ko_no_argument, 'c'},  // apply the logarithmic correction, not fully implemented
                                      {"idle", ko_no_argument, 'i'},
                                      {"help", ko_no_argument, 'h'},
                                      {"prewhiten", ko_no_argument, OPT_PREWHITEN},
                                      {"debug", ko_no_argument, OPT_DEBUG},
                                      // Negative cases, corresponding to long arguments only;
                                      {"strip", ko_no_argument, OPT_STRIP},
                                      {"pack", ko_no_argument, OPT_PACK},
                                      {NULL, 0, 0}};

    // Initialize ketopt, skipping the first two positional arguments
    ketopt_t opt = KETOPT_INIT;
    opt.ind = 2;  // Start parsing options from argv[2]

    int c;
    while ((c = ketopt(&opt, argc, argv, 1, "o:d:n:t:f:e:j:m:g:sich", longopts)) != -1) {
        // printf("argument: %c", c);
        switch (c) {
            case 'o':
                params.oversamplingFactor = atof(opt.arg);
                break;
            case 'n':
                params.npeaks = atoi(opt.arg);
                break;
            case 'd':
                params.nterms = atoi(opt.arg);
                break;
            case 't':
                params.threshold = atof(opt.arg) * M_LN10;
                break;
            case 'f':
                params.fmin = atof(opt.arg);
                break;
            case 'e':
                if (!parse_gb_eval(opt.arg, &params.gbEvalMode, &params.gbAlpha)) {
                    fprintf(stderr, "Invalid evaluation '%s'. Expected gbls, gbl, gbas, or gba with optional [alpha].\n", opt.arg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'j':
                params.jobs = atoi(opt.arg);
                break;
            case 'm':
                errno = 0;
                char *mode_end = NULL;
                long parsed_mode = strtol(opt.arg, &mode_end, 10);
                if (errno != 0 || mode_end == opt.arg || *mode_end != '\0' || parsed_mode < 0 || parsed_mode > 6) {
                    fprintf(stderr, "Invalid mode '%s'. Expected an integer from 0 to 6.\n", opt.arg);
                    exit(EXIT_FAILURE);
                }
                params.mode = (int)parsed_mode;
                break;
            case 'g':
                if (!parse_periodogram_method(opt.arg, &params.periodogramMethod)) {
                    fprintf(stderr, "Invalid periodogram method '%s'. Expected ihs, aov, aovmh, aobmhw, chi, chi2, or fastchi2.\n", opt.arg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 's':
                params.spectrum = true;
                break;
            case OPT_EPSILON:
                params.epsilon = atof(opt.arg);
                break;
            case OPT_NUFFT:
                if (!parse_nufft_mode(opt.arg, &params.gridMode)) {
                    fprintf(stderr, "Invalid NuFFT backend '%s'. Expected 43, pswf43, 21, or pswf21.\n", opt.arg);
                    exit(EXIT_FAILURE);
                }
                break;
            case OPT_PREWHITEN:
                params.prewhiten = true;
                // printf("Prewhitening\n");
                break;
            case OPT_DEBUG:
                params.debug = true;
                break;
            case 'i':
                params.idle = true;
                break;
            case 'c':
                params.corrected = true;
                break;
            case 'h':
                print_help(argv);
                exit(0);
                break;
            case OPT_PACK:
                print_help(argv);
                exit(0);
                break;
            case OPT_STRIP:
                print_help(argv);
                exit(0);
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
static inline double correct_ihs_threshold(const double threshold, const int n) {
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
        if (f_x * f_low < 0)
            high = x;
        else
            low = x;
    }
    return x;
}

static inline double correct_aov_threshold(const double threshold, const int degree, const int n_eff) {
    if (n_eff <= (2 * degree) + 1) return INFINITY;
    double low = 0.0;
    double high = 1.0 - 1.0e-12;
    for (int i = 0; i < 80; ++i) {
        double mid = 0.5 * (low + high);
        double value = lnFAP((2 * degree) + 1, 2 * degree, mid, n_eff);
        if (!double_is_finite_bits(value) || value >= threshold) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return high;
}

static inline double correct_threshold(const parameters *params, const buffer_t *buffer) {
    if (periodogram_uses_aov(params->periodogramMethod)) {
        return correct_aov_threshold(params->threshold, params->nterms, periodogram_effective_n(buffer));
    }
    if (mode_uses_direct_gb_grid(params->mode)) return params->threshold;
    return correct_ihs_threshold(params->threshold, params->nterms);
}

#endif  // PARAMS_H
