#ifndef PROCESS_H
#define PROCESS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readout.h"
#include "params.h"

#include "../include/klib/kstring.h"
#include "../include/mufft/mufft.x86.h"

void process_target(char* in_file, buffer_t* buffer, parameters* params){
    read_dat(kv_A(params->targets, 0).path, buffer); linreg_buffer(buffer); //read the data from .dat file
    int n = 2 + (int)(log(buffer->x[buffer->n-1]) * M_LOG10E); //number of significant digits required for the spectrum

    double number = 0.123456789; printf("\n%i\n", n); // placeholder line
    ksprintf(&buffer->spectrum, "%.*f", n, number); // printf("Formatted string: %s\n", buffer->spectrum.s); //prinft test






    if (params->spectrum){
        // Prepare the output file
        char out_file[256];
        strncpy(out_file, in_file, sizeof(out_file) - 1);
        out_file[sizeof(out_file) - 1] = '\0'; // Ensure null-termination
        char* dot = strrchr(out_file, '.'); // Find the last dot in the file name
        if (dot) {*dot = '\0';}  // Remove the existing extension
        strcat(out_file, ".tsv"); // Append the new extension

        // Write the formatted string to the new .tsv file
        FILE* fp = fopen(out_file, "w");
        if (fp == NULL) {perror("Failed to open file for writing"); return;}

        fprintf(fp, "%s\n", buffer->spectrum.s); fclose(fp); // Write the spectrum string to the file
    }

    if(buffer->spectrum.s){free(buffer->spectrum.s);}
}



#endif
