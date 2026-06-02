#ifndef STRINGS_H
#define STRINGS_H

#include <fast_convert.h>
#include <sds.h>

inline unsigned int custom_ftoa(float v, int size, char *line) {
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1) {
        line[0] = '0';
        line[1] = '.';
        for (int i = 0; i < size; i++) {
            line[2 + i] = '0';
        }
        line[4] = '\0';
        return 4;  // Return the length of the string
    } else {
        return fast_ftoa(v, new_size, line);
    }
}

inline unsigned int custom_dtoa(double v, int size, char *line) {
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1) {
        return fast_dtoa(0.0, new_size, line);
    } else {
        return fast_dtoa(v, new_size, line);
    }
}

static inline double output_coordinate(double freq, bool outputPeriod) {
    if (!outputPeriod) return freq;
    if (freq > 0.0 && double_is_finite_bits(freq)) return 1.0 / freq;
    return INFINITY;
}

static inline int output_frequency_significant_digits(double freq, int precision) {
    if (!(freq > 0.0) || !double_is_finite_bits(freq)) return precision;
    int digits = precision + (int)ceil(log10(freq));
    if (digits < 1) digits = 1;
    if (digits > PREC_DBL_NR) digits = PREC_DBL_NR;
    return digits;
}

static inline int output_fractional_digits(double coord, int significant_digits) {
    if (!(coord > 0.0) || !double_is_finite_bits(coord)) return 0;
    int integer_digits = (int)ceil(log10(coord));
    if (integer_digits < 0) integer_digits = 0;
    int fractional_digits = significant_digits - integer_digits;
    return fractional_digits > 0 ? fractional_digits : 0;
}

static inline void pad_fractional_digits(char *line, int fractional_digits) {
    if (fractional_digits <= 0 || strchr(line, 'e') || strchr(line, 'E')) return;
    char *dot = strchr(line, '.');
    size_t len = strlen(line);
    if (!dot) {
        line[len++] = '.';
        line[len] = '\0';
        dot = line + len - 1U;
    }
    int current = (int)(len - (size_t)(dot - line) - 1U);
    while (current < fractional_digits) {
        line[len++] = '0';
        line[len] = '\0';
        ++current;
    }
}

static inline void custom_coordinate_dtoa(double freq, bool outputPeriod, int precision, char *line) {
    double coord = output_coordinate(freq, outputPeriod);
    if (!double_is_finite_bits(coord)) {
        strcpy(line, "inf");
        return;
    }
    if (outputPeriod) {
        int significant_digits = output_frequency_significant_digits(freq, precision);
        fast_dtoa(coord, significant_digits, line);
        pad_fractional_digits(line, output_fractional_digits(coord, significant_digits));
        return;
    }
    custom_dtoa(coord, precision, line);
}

static inline void appendFreq(double freq, float magnitude, int precision, bool outputPeriod, sds *buffer, char *stringBuff) {
    // Convert output coordinate to string and append to the buffer.
    custom_coordinate_dtoa(freq, outputPeriod, precision, stringBuff);
    *buffer = sdscat(*buffer, stringBuff);
    *buffer = sdscatlen(*buffer, "\t", 1);  // Add the separator

    // Convert magnitude to string and append to the buffer
    custom_ftoa(magnitude, 2, stringBuff);
    *buffer = sdscat(*buffer, stringBuff);
    *buffer = sdscatlen(*buffer, "\n", 1);  // Add a newline character
}

static inline void write_tsv(buffer_t *buffer, char *in_file) {
    char out_file[256];
    strncpy(out_file, in_file, sizeof(out_file) - 1);
    out_file[sizeof(out_file) - 1] = '\0';  // Ensure null-termination
    char *dot = strrchr(out_file, '.');     // Find the last dot in the file name
    if (dot) {
        *dot = '\0';
    }  // Remove the existing extension
    strcat(out_file, ".tsv");  // Append the new extension

    // Write the formatted string to the new .tsv file
    FILE *fp = fopen(out_file, "w");
    if (fp == NULL) {
        perror("Failed to open file for writing");
        return;
    }

    fprintf(fp, "%s\n", buffer->spectrum);
    fclose(fp);  // Write the spectrum string to the file
    sdsclear(buffer->spectrum);
};

static inline void write_spectrum_matrix_tsv(char *in_file, double fmin, double fstep, uint32_t nfreq, uint32_t ncols, const float *matrix, int precision,
                                             bool outputPeriod, char *stringBuff) {
    if (!matrix || nfreq == 0U || ncols == 0U) return;

    char out_file[256];
    strncpy(out_file, in_file, sizeof(out_file) - 1);
    out_file[sizeof(out_file) - 1] = '\0';
    char *dot = strrchr(out_file, '.');
    if (dot) {
        *dot = '\0';
    }
    strcat(out_file, ".tsv");

    FILE *fp = fopen(out_file, "w");
    if (fp == NULL) {
        perror("Failed to open file for writing");
        return;
    }

    for (uint32_t i = 0; i < nfreq; ++i) {
        double freq = fmin + ((double)i * fstep);
        custom_coordinate_dtoa(freq, outputPeriod, precision, stringBuff);
        fputs(stringBuff, fp);
        for (uint32_t col = 0; col < ncols; ++col) {
            fputc('\t', fp);
            custom_ftoa(matrix[((size_t)col * (size_t)nfreq) + (size_t)i], 2, stringBuff);
            fputs(stringBuff, fp);
        }
        fputc('\n', fp);
    }

    fclose(fp);
}

static inline void append_spaces(sds *buffer, int count) {
    char spaces[64];
    while (count > 0) {
        int chunk = count > (int)sizeof(spaces) ? (int)sizeof(spaces) : count;
        memset(spaces, ' ', (size_t)chunk);
        *buffer = sdscatlen(*buffer, spaces, (size_t)chunk);
        count -= chunk;
    }
}

static inline void append_center_cell(sds *buffer, const char *str, int width) {
    int len = (int)strlen(str);
    int pad = width - len;
    int left = pad > 0 ? pad / 2 : 0;
    int right = pad > 0 ? pad - left : 0;
    append_spaces(buffer, left);
    *buffer = sdscat(*buffer, str);
    append_spaces(buffer, right);
}

void print_peaks(buffer_t *buffer, parameters *params, int n, char *stringBuff, char *in_file, int mode, gb_eval_mode evalMode) {
    if (mode > 1) {
        n += 1;
    }
    bool use_aov = periodogram_uses_aov(params->periodogramMethod);
    bool print_amp = mode > 0;
    const eval_method_t *method = eval_method_for_mode(evalMode);
    bool print_width = print_amp && method->width_label;
    // Column indices: 0: output coordinate, 1: NLL, 2: Amp/R2, 3: statistic, 4: BLS duration.
    const char *hdr[5] = {params->outputPeriod ? "P[d]" : "f[1/d]", method->peak_power_label, "Amp", method->stat_label, method->width_label};
    if (use_aov && !print_amp) hdr[2] = "R2";
    int HdrCols = print_amp ? (print_width ? 5 : 4) : (use_aov ? 3 : 2);
    int colWidth[5];
    int i;

    for (i = 0; i < HdrCols; i++) colWidth[i] = (int)strlen(hdr[i]);
    for (i = 0; i < buffer->nPeaks && buffer->peaks[i].p > 0; i++) {
        custom_coordinate_dtoa(buffer->peaks[i].freq, params->outputPeriod, n, stringBuff);
        int len = (int)strlen(stringBuff);
        if (len > colWidth[0]) colWidth[0] = len;

        custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
        len = (int)strlen(stringBuff);
        if (len > colWidth[1]) colWidth[1] = len;

        if (use_aov && !print_amp) {
            custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
            len = (int)strlen(stringBuff);
            if (len > colWidth[2]) colWidth[2] = len;
        } else if (print_amp) {
            custom_ftoa(buffer->peaks[i].amp, 3, stringBuff);
            len = (int)strlen(stringBuff);
            if (len > colWidth[2]) colWidth[2] = len;

            custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
            len = (int)strlen(stringBuff);
            if (len > colWidth[3]) colWidth[3] = len;

            if (print_width) {
                custom_ftoa(buffer->peaks[i].width, 2, stringBuff);
                len = (int)strlen(stringBuff);
                if (len > colWidth[4]) colWidth[4] = len;
            }
        }
    }

    buffer->outBuf = sdscat(buffer->outBuf, "ihsnpeaks ");
    buffer->outBuf = sdscat(buffer->outBuf, IHSNPEAKS_VERSION);
    buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);

    // Append file information (unchanged)
    buffer->outBuf = sdscat(buffer->outBuf, "File: ");
    buffer->outBuf = sdscat(buffer->outBuf, in_file);
    buffer->outBuf = sdscat(buffer->outBuf, "    NofData: ");
    custom_ftoa(buffer->n, 0, stringBuff);  // using custom_ftoa to print integer
    buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
    buffer->outBuf = sdscat(buffer->outBuf, "    Average: ");
    custom_ftoa(buffer->magnitude, 2, stringBuff);
    buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
    buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);

    // Print table header.
    {
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        for (i = 0; i < HdrCols; i++) {
            append_center_cell(&buffer->outBuf, hdr[i], colWidth[i]);
            if (i + 1 < HdrCols)
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
            else
                buffer->outBuf = sdscat(buffer->outBuf, "\n");
        }
    }
    // Second pass: print each peak row with proper padding.
    {
        for (i = 0; i < buffer->nPeaks && buffer->peaks[i].p > 0; i++) {
            buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);

            custom_coordinate_dtoa(buffer->peaks[i].freq, params->outputPeriod, n, stringBuff);
            append_center_cell(&buffer->outBuf, stringBuff, colWidth[0]);
            buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);

            custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
            append_center_cell(&buffer->outBuf, stringBuff, colWidth[1]);

            if (use_aov && !print_amp) {
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
                custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
                append_center_cell(&buffer->outBuf, stringBuff, colWidth[2]);
            } else if (print_amp) {
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
                custom_ftoa(buffer->peaks[i].amp, 3, stringBuff);
                append_center_cell(&buffer->outBuf, stringBuff, colWidth[2]);

                buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
                custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
                append_center_cell(&buffer->outBuf, stringBuff, colWidth[3]);

                if (print_width) {
                    buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
                    custom_ftoa(buffer->peaks[i].width, 2, stringBuff);
                    append_center_cell(&buffer->outBuf, stringBuff, colWidth[4]);
                }
            }
            buffer->outBuf = sdscat(buffer->outBuf, "\n");
        }
    }

    printf("%s", buffer->outBuf);
    sdsclear(buffer->outBuf);
}

void fprint_buffer(buffer_t *buffer, parameters *params) {
    pthread_mutex_lock(&params->mutex);
    FILE *file = fopen(params->outFile, "a");  // Open the file in append mode
    if (file == NULL) {
        pthread_mutex_unlock(&params->mutex);
        perror("Failed to open file for appending");
        return;
    }
    fprintf(file, "%s", buffer->outBuf);
    fflush(stdout);
    fclose(file);
    sdsclear(buffer->outBuf);
    pthread_mutex_unlock(&params->mutex);
}

void append_peaks(buffer_t *buffer, parameters *params, int n, char *stringBuff, char *in_file, int mode, gb_eval_mode evalMode) {
    const eval_method_t *method = eval_method_for_mode(evalMode);
    bool print_width = mode > 0 && method->width_label;
    bool use_aov = periodogram_uses_aov(params->periodogramMethod);
    bool print_amp = mode > 0;
    int i = 0;
    if (buffer->nPeaks > 0) {
        // Append file information to the output buffer
        buffer->outBuf = sdscat(buffer->outBuf, in_file);
        custom_coordinate_dtoa(buffer->peaks[0].freq, params->outputPeriod, n, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        custom_ftoa(buffer->peaks[0].p * M_LOG10E, 2, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        // buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);

        // Append each peak's information to the output buffer
        if (use_aov && !print_amp) {
            while (i < buffer->nPeaks && buffer->peaks[i].p > 0) {
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
                custom_coordinate_dtoa(buffer->peaks[i].freq, params->outputPeriod, n + (mode > 1 ? 1 : 0), stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
                i++;
            }
        } else if (use_aov) {
            while (i < buffer->nPeaks && buffer->peaks[i].p > 0) {
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
                custom_coordinate_dtoa(buffer->peaks[i].freq, params->outputPeriod, n + 1, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                custom_ftoa(buffer->peaks[i].amp, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                if (print_width) {
                    buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                    custom_ftoa(buffer->peaks[i].width, 2, stringBuff);
                    buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                }
                buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
                i++;
            }
        } else if (mode == 0) {
            while (i < buffer->nPeaks && buffer->peaks[i].p > 0) {
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
                custom_coordinate_dtoa(buffer->peaks[i].freq, params->outputPeriod, n, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                // Convert base-10 NLL to string using custom_ftoa
                custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
                i++;
            }
        } else {
            while (i < buffer->nPeaks && buffer->peaks[i].p > 0) {
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
                custom_coordinate_dtoa(buffer->peaks[i].freq, params->outputPeriod, n + 1, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                // Convert amplitude to string using custom_ftoa
                custom_ftoa(buffer->peaks[i].amp, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                custom_ftoa(buffer->peaks[i].r2, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                if (print_width) {
                    buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                    custom_ftoa(buffer->peaks[i].width, 2, stringBuff);
                    buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                }
                buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
                i++;
            }
        }
        buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);
        buffer->loc_iter += 1;
        if (buffer->loc_iter == 0) {  // append to the file only if the local counter oveflows
            // Append the buffer content to the file
            fprint_buffer(buffer, params);
            fflush(stdout);
        }
    }
}

#endif
