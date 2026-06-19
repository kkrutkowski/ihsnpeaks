#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include "../../src/nufft/nufft1.h"
#include "../src/nufft/cunufft.h"

double get_time_ms() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec * 1e-6;
}

void generate_data(int M, double *x, float *y_real, float *y_imag) {
    for (int m = 0; m < M; ++m) {
        x[m] = (double)rand() / RAND_MAX;
        y_real[m] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
        y_imag[m] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
}

int get_physical_cores() {
    FILE *fp = fopen("/proc/cpuinfo", "r");
    if (!fp) return 1;
    char line[256];
    int cores = 1;
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "cpu cores", 9) == 0) {
            char *colon = strchr(line, ':');
            if (colon) {
                cores = atoi(colon + 1);
                break;
            }
        }
    }
    fclose(fp);
    return cores > 0 ? cores : 1;
}

void run_benchmark(cunufft_mode mode_type) {
    const int Mpoints = 1000;
    const int freq_factor = 32;
    const double df = 1.0;

    nufft1_mode cpu_mode = (mode_type == CUNUFFT_PSWF21) ? NUFFT1_PSWF21 : NUFFT1_PSWF43;
    cunufft_mode gpu_mode = mode_type;

    /* ---- Warm up the OpenMP thread pool (avoid first-call overhead) ---- */
    #pragma omp parallel
    { volatile int dummy = omp_get_thread_num(); (void)dummy; }

    int physical_cores = get_physical_cores();
    int num_threads = 1; // Run CPU on a single thread to avoid memory bus contention

    printf("\n=====================================================================================================================\n");
    printf("              M = %d points, %d transform batches, Mode = %s (CPU Est. scaled by %d physical cores)\n",
           Mpoints, freq_factor,
           (mode_type == CUNUFFT_PSWF21) ? "PSWF21" : "PSWF43",
           physical_cores);
    printf("=====================================================================================================================\n");
    printf("%-7s | %-8s | %-14s | %-14s | %-14s | %-10s | %-10s | %-10s\n",
           "Log2(N)", "Nfft", "CPU Est (ms)", "GPU Kern (ms)", "GPU Total (ms)", "Speedup(D)", "Max Err", "Avg Err");
    printf("---------------------------------------------------------------------------------------------------------------------\n");

    double *x = (double *)malloc(Mpoints * sizeof(double));
    float *y_real = (float *)malloc(Mpoints * sizeof(float));
    float *y_imag = (float *)malloc(Mpoints * sizeof(float));

    generate_data(Mpoints, x, y_real, y_imag);

    for (int log2N = 10; log2N <= 20; ++log2N) {
        int N = 1 << log2N;
        if (cpu_mode == NUFFT1_PSWF43 && (N % 4) != 0) continue;

        /* =========== CPU NuFFT Setup =========== */
        nufft1_external_sizes cpu_sizes = nufft1_external_sizes_for_plan(N, cpu_mode);

        float *cpu_tw_r = (float *)aligned_alloc(64, cpu_sizes.twiddle_len * sizeof(float));
        float *cpu_tw_i = (float *)aligned_alloc(64, cpu_sizes.twiddle_len * sizeof(float));

        nufft1_plan *cpu_plan = nufft1_initialize_shared(Mpoints, N, freq_factor, cpu_mode, cpu_tw_r, cpu_tw_i);
        if (!cpu_plan) {
            printf("Failed to initialize CPU plan for N = %d\n", N);
            free(cpu_tw_r); free(cpu_tw_i);
            continue;
        }

        /*
         * Create one workspace PER THREAD (not per factor).
         * Each workspace stores precomputed spread data for all 32 factors
         * but has its own private FFT/cobra buffers for thread-safe execute.
         * This reduces the memory footprint from 32 * 18 MB = 576 MB
         * down to num_threads * 18 MB, and critically avoids 32 concurrent
         * FFT buffers thrashing the L3 cache.
         */
        nufft1_workspace **tws    = (nufft1_workspace **)malloc(num_threads * sizeof(*tws));
        float            **t_fr   = (float **)malloc(num_threads * sizeof(float*));
        float            **t_fi   = (float **)malloc(num_threads * sizeof(float*));
        float            **t_cr   = (float **)malloc(num_threads * sizeof(float*));
        float            **t_ci   = (float **)malloc(num_threads * sizeof(float*));

        for (int t = 0; t < num_threads; ++t) {
            t_fr[t] = (float *)aligned_alloc(64, cpu_sizes.fft_len   * sizeof(float));
            t_fi[t] = (float *)aligned_alloc(64, cpu_sizes.fft_len   * sizeof(float));
            t_cr[t] = (float *)aligned_alloc(64, cpu_sizes.cobra_len * sizeof(float));
            t_ci[t] = (float *)aligned_alloc(64, cpu_sizes.cobra_len * sizeof(float));

            /*
             * Touch all pages NOW so that the kernel maps physical pages before
             * timing. Without this, the first nufft1_execute() memset triggers
             * page faults, and 32 threads page-faulting concurrently serialise
             * on the kernel's page-fault handler lock.
             */
            memset(t_fr[t], 0, cpu_sizes.fft_len   * sizeof(float));
            memset(t_fi[t], 0, cpu_sizes.fft_len   * sizeof(float));
            memset(t_cr[t], 0, cpu_sizes.cobra_len * sizeof(float));
            memset(t_ci[t], 0, cpu_sizes.cobra_len * sizeof(float));

            tws[t] = nufft1_create_workspace(cpu_plan, t_fr[t], t_fi[t], t_cr[t], t_ci[t]);
            nufft1_precompute(tws[t], x, Mpoints, df);
        }

        /* =========== GPU NuFFT Setup =========== */
        cunufft_plan *gpu_plan = cunufft_initialize(Mpoints, N, freq_factor, gpu_mode);
        if (!gpu_plan) {
            printf("Failed to initialize GPU plan for N = %d\n", N);
            for (int t = 0; t < num_threads; ++t) {
                nufft1_free_workspace(tws[t]);
                free(t_fr[t]); free(t_fi[t]);
                free(t_cr[t]); free(t_ci[t]);
            }
            free(tws); free(t_fr); free(t_fi); free(t_cr); free(t_ci);
            nufft1_free_plan(cpu_plan);
            free(cpu_tw_r); free(cpu_tw_i);
            continue;
        }
        cunufft_workspace *gpu_ws = cunufft_create_workspace(gpu_plan);
        cunufft_precompute(gpu_ws, x, Mpoints, df);

        /* =========== Output buffers =========== */
        int Nout = nufft1_output_len(cpu_plan);
        size_t total_out_bytes = (size_t)freq_factor * (size_t)Nout * sizeof(float);
        float *cpu_out_r = (float *)aligned_alloc(64, total_out_bytes);
        float *cpu_out_i = (float *)aligned_alloc(64, total_out_bytes);
        float *gpu_out_r = (float *)aligned_alloc(64, total_out_bytes);
        float *gpu_out_i = (float *)aligned_alloc(64, total_out_bytes);

        /* Pre-fault output pages */
        memset(cpu_out_r, 0, total_out_bytes);
        memset(cpu_out_i, 0, total_out_bytes);

        /* =========== CPU warm-up & timed run (single parallel block to avoid thread wakeup overhead) =========== */
        double cpu_start = 0.0, cpu_end = 0.0;
        #pragma omp parallel num_threads(num_threads)
        {
            int tid = omp_get_thread_num();

            /* 1. Warm-up run */
            #pragma omp for schedule(static)
            for (int f = 1; f <= freq_factor; ++f)
                nufft1_execute(tws[tid], y_real, y_imag,
                               cpu_out_r + (f-1)*Nout,
                               cpu_out_i + (f-1)*Nout, f);

                #pragma omp barrier
                #pragma omp master
                {
                    cpu_start = get_time_ms();
                }
                #pragma omp barrier

                /* 2. Timed run */
                #pragma omp for schedule(static)
                for (int f = 1; f <= freq_factor; ++f)
                    nufft1_execute(tws[tid], y_real, y_imag,
                                   cpu_out_r + (f-1)*Nout,
                                   cpu_out_i + (f-1)*Nout, f);

                    #pragma omp barrier
                    #pragma omp master
                    {
                        cpu_end = get_time_ms();
                    }
        }
        double cpu_time = (cpu_end - cpu_start) / (double)physical_cores;

        /* =========== GPU warm-up & timed run =========== */
        cunufft_execute_batch(gpu_ws, y_real, y_imag, gpu_out_r, gpu_out_i);

        double gpu_total_start = get_time_ms();
        float gpu_kernel_ms = cunufft_execute_batch(gpu_ws, y_real, y_imag, gpu_out_r, gpu_out_i);
        double gpu_total_end = get_time_ms();
        double gpu_total_time = gpu_total_end - gpu_total_start;

        /* =========== Error Calculation =========== */
        double max_err = 0.0, avg_err = 0.0;
        int count = 0;
        for (int f = 1; f <= freq_factor; ++f) {
            int offset = (f - 1) * Nout;
            for (int k = 0; k < Nout; ++k) {
                double dr = (double)cpu_out_r[offset+k] - (double)gpu_out_r[offset+k];
                double di = (double)cpu_out_i[offset+k] - (double)gpu_out_i[offset+k];
                double err = sqrt(dr*dr + di*di);
                if (err > max_err) max_err = err;
                avg_err += err;
                ++count;
            }
        }
        avg_err /= count;

        double speedup_device = cpu_time / (double)gpu_kernel_ms;

        printf("2^%-5d | %-8d | %-14.4f | %-14.4f | %-14.4f | %-10.2fx | %-10.4e | %-10.4e\n",
               log2N, (int)cpu_sizes.twiddle_len,
               cpu_time, (double)gpu_kernel_ms, gpu_total_time,
               speedup_device, max_err, avg_err);

        /* =========== Cleanup =========== */
        free(cpu_out_r); free(cpu_out_i);
        free(gpu_out_r); free(gpu_out_i);

        for (int t = 0; t < num_threads; ++t) {
            nufft1_free_workspace(tws[t]);
            free(t_fr[t]); free(t_fi[t]);
            free(t_cr[t]); free(t_ci[t]);
        }
        free(tws); free(t_fr); free(t_fi); free(t_cr); free(t_ci);

        nufft1_free_plan(cpu_plan);
        free(cpu_tw_r); free(cpu_tw_i);

        cunufft_free_workspace(gpu_ws);
        cunufft_free_plan(gpu_plan);
    }
    printf("=====================================================================================================================\n");
    free(x); free(y_real); free(y_imag);
}

int main() {
    setenv("OMP_WAIT_POLICY", "passive", 1);
    srand(42);
    run_benchmark(CUNUFFT_PSWF21);
    run_benchmark(CUNUFFT_PSWF43);
    return 0;
}
