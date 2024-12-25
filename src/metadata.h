#ifndef METADATA_H
#define METADATA_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <stdbool.h>

#include "params.h"
#include "../include/klib/kvec.h"

// Define the target struct to hold file path and a pointer to the parameters struct
typedef struct {
    char *path;
    parameters *params;
} target_t;

// Define a kvec type for target_t structs
typedef kvec_t(target_t) kvec_target_t;

// Function to check if the path is a file or directory
static int is_directory(const char *path) {
    struct stat path_stat;
    if (stat(path, &path_stat) != 0) {
        perror("stat");
        exit(EXIT_FAILURE);
    }
    return S_ISDIR(path_stat.st_mode);
}

// Function to process the path and populate the target vector
static void process_path(parameters *params, kvec_target_t *targets, uint32_t *maxLen) {
    kv_init(*targets);
    *maxLen = 0;

    const char *path = params->target[0]; // Use the first target path

    if (!is_directory(path)) {
        params->isFile = true;

        // Allocate and populate a target struct with the file path and params pointer
        target_t target;
        target.path = strdup(path);
        target.params = params;
        kv_push(target_t, *targets, target); // Store the target in the vector

        // Check newline count for this single file
        FILE *file = fopen(path, "r");
        if (!file) {
            perror("fopen");
            exit(EXIT_FAILURE);
        }

        uint32_t newline_count = 0;
        int ch;
        while ((ch = fgetc(file)) != EOF) {
            if (ch == '\n') {
                newline_count++;
            }
        }
        fclose(file);

        *maxLen = newline_count;
    } else {
        params->isFile = false;
        DIR *dir = opendir(path);
        if (!dir) {
            perror("opendir");
            exit(EXIT_FAILURE);
        }

        struct dirent *entry;
        off_t largest_size = 0;

        while ((entry = readdir(dir)) != NULL) {
            if (entry->d_type == DT_REG) {
                // Construct the full path for each file
                char full_path[256];
                snprintf(full_path, sizeof(full_path), "%s/%s", path, entry->d_name);

                // Get the file size
                struct stat file_stat;
                if (stat(full_path, &file_stat) != 0) {
                    perror("stat");
                    exit(EXIT_FAILURE);
                }

                // Check if it's the largest file so far
                if (file_stat.st_size > largest_size) {
                    largest_size = file_stat.st_size;

                    // Open and count newlines in this file
                    FILE *file = fopen(full_path, "r");
                    if (!file) {
                        perror("fopen");
                        exit(EXIT_FAILURE);
                    }

                    uint32_t newline_count = 0;
                    int ch;
                    while ((ch = fgetc(file)) != EOF) {
                        if (ch == '\n') {
                            newline_count++;
                        }
                    }
                    fclose(file);

                    *maxLen = newline_count;
                }

                // Allocate and populate a target struct with the file path and params pointer
                target_t target;
                target.path = strdup(full_path);
                target.params = params;
                kv_push(target_t, *targets, target); // Store the target in the vector
            }
        }
        closedir(dir);
    }
}

// Function to free all memory allocated in the kvec and destroy it
static void free_targets(kvec_target_t *targets) {
    for (size_t i = 0; i < kv_size(*targets); i++) {
        free(kv_A(*targets, i).path);
    }
    kv_destroy(*targets);
}

#endif // METADATA_H
