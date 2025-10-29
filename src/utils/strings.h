#ifndef STRINGS_H
#define STRINGS_H

#include <sds.h>
#include <fast_convert.h>


inline unsigned int custom_ftoa(float v, int size, char *line) {
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1) {
        line[0] = '0'; line[1] = '.';
        for (int i = 0; i < size; i++) {line[2 + i] = '0';}
        line[4] = '\0'; return 4; // Return the length of the string
    } else {return fast_ftoa(v, new_size, line);}
}

inline unsigned int custom_dtoa(double v, int size, char *line){
    int new_size = size + (int)ceil(log10(v));
    if (new_size < 1){return fast_dtoa(0.0, new_size, line);}
    else {return fast_dtoa(v, new_size, line);}
}

static inline void appendFreq(double freq, float magnitude, int precision, sds* buffer, char* stringBuff) {
    // Convert frequency to string and append to the buffer
    custom_dtoa(freq, precision, stringBuff);
    *buffer = sdscat(*buffer, stringBuff);
    *buffer = sdscatlen(*buffer, "\t", 1);  // Add the separator

    // Convert magnitude to string and append to the buffer
    custom_ftoa(magnitude, 2, stringBuff);
    *buffer = sdscat(*buffer, stringBuff);
    *buffer = sdscatlen(*buffer, "\n", 1);  // Add a newline character
}

static inline void write_tsv(buffer_t *buffer, char* in_file){
        char out_file[256];
        strncpy(out_file, in_file, sizeof(out_file) - 1);
        out_file[sizeof(out_file) - 1] = '\0'; // Ensure null-termination
        char* dot = strrchr(out_file, '.'); // Find the last dot in the file name
        if (dot) {*dot = '\0';}  // Remove the existing extension
        strcat(out_file, ".tsv"); // Append the new extension

        // Write the formatted string to the new .tsv file
        FILE* fp = fopen(out_file, "w");
        if (fp == NULL) {perror("Failed to open file for writing"); return;}

        fprintf(fp, "%s\n", buffer->spectrum); fclose(fp); // Write the spectrum string to the file
        sdsclear(buffer->spectrum);
};

static int column_width(double val, int precision) {
    int intPartWidth = (val > 0) ? (int)ceil(log10(val)) : 1;
    return 3 + precision + intPartWidth;
}

void print_peaks(buffer_t *buffer, parameters *params, int n, char *stringBuff, char *in_file, int mode) {
    if (mode > 1){n += 1;}
    // Column indices: 0: f[1/d] (frequency), 1: log(p), 2: Amp, 3: R.
    const char *hdr[4] = { "f[1/d]", "log(p)", "Amp", "R" };
    int hdrWidth[4]; int colWidth[4]; int i;

    // Set header widths (without any extra padding)
    for(i = 0; i < 4; i++){hdrWidth[i] = (int)strlen(hdr[i]); colWidth[i] = hdrWidth[i];}

    /* First pass: scan through peaks and compute the maximum printed width for each column.
       Note: For log(p) we compute log10(p) by first converting the natural log to log10 using M_LOG10E. */
    for(i = 0; i < buffer->nPeaks && buffer->peaks[i].p > 0; i++){
        double freq = buffer->peaks[i].freq;
        double logp = buffer->peaks[i].p * M_LOG10E; // equivalent to log10(p)
        double amp  = buffer->peaks[i].amp;
        double r    = buffer->peaks[i].r;
        //printf("%.3f\n", r); // ??? - wrong here
        //if (r <= 1.0 && mode > 0) {r = get_r(buffer, buffer->peaks[i].freq, NULL, false);}
        int w;
        w = column_width(freq, n); if (w > colWidth[0]) colWidth[0] = w;
        w = column_width(logp, 2); if (w > colWidth[1]) colWidth[1] = w;
        w = column_width(amp, 3); if (w > colWidth[2]) colWidth[2] = w;
        w = column_width(r, 2); if (w > colWidth[3]) colWidth[3] = w;
    }

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

    // Print table header with each header padded on the left to the column width.
    {
        char temp[64]; // temporary buffer for spaces
        int pad; int HdrCols = 4; if (mode == 0){HdrCols = 2;}
        for (i = 0; i < HdrCols; i++) {
            pad = colWidth[i] - (int)strlen(hdr[i]);
            if (pad > 0) {
                memset(temp, ' ', pad);
                temp[pad] = '\0';
                buffer->outBuf = sdscat(buffer->outBuf, temp);
            }
            buffer->outBuf = sdscat(buffer->outBuf, hdr[i]);
            if (i < 3) buffer->outBuf = sdscat(buffer->outBuf, "\t");
            else buffer->outBuf = sdscat(buffer->outBuf, "\n");
        }
    }
    if (mode == 0){buffer->outBuf = sdscat(buffer->outBuf, "\n");}
    // Second pass: print each peak row with proper padding.
    {
        char temp[64];
        int pad, len;
        for (i = 0; i < buffer->nPeaks && buffer->peaks[i].p > 0; i++){
            // Frequency column (precision n)
            custom_dtoa(buffer->peaks[i].freq, n, stringBuff);
            len = (int)strlen(stringBuff);
            pad = colWidth[0] - len;
            if (pad > 0) {
                memset(temp, ' ', pad);
                temp[pad] = '\0';
                buffer->outBuf = sdscat(buffer->outBuf, temp);
            }
            buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
            buffer->outBuf = sdscat(buffer->outBuf, "\t");

            // log(p) (precision 2)
            custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
            len = (int)strlen(stringBuff);
            pad = colWidth[1] - len;
            if (pad > 0) {
                memset(temp, ' ', pad);
                temp[pad] = '\0';
                buffer->outBuf = sdscat(buffer->outBuf, temp);
            }
            buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
            buffer->outBuf = sdscat(buffer->outBuf, "\t");

            if (mode > 0){
                // Amplitude (precision 3)
                custom_ftoa(buffer->peaks[i].amp, 3, stringBuff);
                len = (int)strlen(stringBuff);
                pad = colWidth[2] - len;
                if (pad > 0) {
                    memset(temp, ' ', pad);
                    temp[pad] = '\0';
                    buffer->outBuf = sdscat(buffer->outBuf, temp);
                }
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, "\t");

                // Test statistics (precision 3)
                custom_ftoa(buffer->peaks[i].r, 3, stringBuff);
                len = (int)strlen(stringBuff);
                pad = colWidth[3] - len;
                if (pad > 0) {
                    memset(temp, ' ', pad);
                    temp[pad] = '\0';
                    buffer->outBuf = sdscat(buffer->outBuf, temp);
                }
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
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
    if (file == NULL) {pthread_mutex_unlock(&params->mutex); perror("Failed to open file for appending"); return;}
    fprintf(file, "%s", buffer->outBuf); fflush(stdout); fclose(file);
    sdsclear(buffer->outBuf);
    pthread_mutex_unlock(&params->mutex);
}

void append_peaks(buffer_t *buffer, parameters *params, int n, char *stringBuff, char *in_file, int mode) {
    int i = 0;
        if (buffer->nPeaks > 0){
        // Append file information to the output buffer
        buffer->outBuf = sdscat(buffer->outBuf, in_file);
        custom_dtoa(buffer->peaks[0].freq, n, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        custom_ftoa(buffer->peaks[0].p * M_LOG10E, 2, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
        buffer->outBuf = sdscatlen(buffer->outBuf, "\t", 1);
        //buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);

        // Append each peak's information to the output buffer
        if (mode == 0){
            while (i < buffer->nPeaks && buffer->peaks[i].p > 0) {
                // Convert frequency to string using custom_ftoa
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
                custom_dtoa(buffer->peaks[i].freq, n, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                // Convert log(p) to string using custom_ftoa
                custom_ftoa(buffer->peaks[i].p * M_LOG10E, 2, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
                i++;
        }} else {
            while (i < buffer->nPeaks && buffer->peaks[i].p > 0) {
                // Convert frequency to string using custom_ftoa
                buffer->outBuf = sdscatlen(buffer->outBuf, "\t[", 2);
                custom_dtoa(buffer->peaks[i].freq, n+1, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                // Convert amplitude to string using custom_ftoa
                custom_ftoa(buffer->peaks[i].amp, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, ", ", 2);
                custom_ftoa(buffer->peaks[i].r, 3, stringBuff);
                buffer->outBuf = sdscat(buffer->outBuf, stringBuff);
                buffer->outBuf = sdscatlen(buffer->outBuf, "]", 1);
                i++;
        }}
        buffer->outBuf = sdscatlen(buffer->outBuf, "\n", 1);
        buffer->loc_iter += 1;
        if(buffer->loc_iter == 0){ //append to the file only if the local counter oveflows
            // Append the buffer content to the file
            fprint_buffer(buffer, params); fflush(stdout);
        }
    }
}

#endif
