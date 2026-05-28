#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "nufft1.h"

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static void fill_inputs(double *x, float *yr, float *yi, int m) {
    srand(1);
    for (int i = 0; i < m; ++i) {
        x[i] = (double)rand() / (double)RAND_MAX;
        yr[i] = (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
        yi[i] = (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
    }
}

static size_t round64(size_t size) { return (size + 63U) & ~(size_t)63U; }

static float *alloc_aligned_float_local(size_t n) {
    float *ptr = (float *)aligned_alloc(64, round64(n * sizeof(float)));
    if (ptr) memset(ptr, 0, round64(n * sizeof(float)));
    return ptr;
}

static double time_precompute(nufft1_mode mode, const double *x, int m, int n, int reps) {
    double total = 0.0;
    for (int r = 0; r < reps; ++r) {
        nufft1_external_sizes sizes = nufft1_external_sizes_for_plan(n, mode);
        float *tw_r = alloc_aligned_float_local(sizes.twiddle_len);
        float *tw_i = alloc_aligned_float_local(sizes.twiddle_len);
        float *fft_r = alloc_aligned_float_local(sizes.fft_len);
        float *fft_i = alloc_aligned_float_local(sizes.fft_len);
        float *cobra_r = alloc_aligned_float_local(sizes.cobra_len);
        float *cobra_i = alloc_aligned_float_local(sizes.cobra_len);
        nufft1_plan *plan = nufft1_initialize_shared(m, n, 1, mode, tw_r, tw_i);
        nufft1_workspace *workspace = nufft1_create_workspace(plan, fft_r, fft_i, cobra_r, cobra_i);
        if (!plan || !workspace) return NAN;
        double t0 = now_seconds();
        nufft1_precompute(workspace, x, m, 1.0);
        total += now_seconds() - t0;
        nufft1_free_workspace(workspace);
        nufft1_free_plan(plan);
        free(tw_r);
        free(tw_i);
        free(fft_r);
        free(fft_i);
        free(cobra_r);
        free(cobra_i);
    }
    return total / (double)reps;
}

static double time_execute(nufft1_mode mode, const double *x, const float *yr, const float *yi, int m, int n, int nout, int reps) {
    nufft1_external_sizes sizes = nufft1_external_sizes_for_plan(n, mode);
    float *tw_r = alloc_aligned_float_local(sizes.twiddle_len);
    float *tw_i = alloc_aligned_float_local(sizes.twiddle_len);
    float *fft_r = alloc_aligned_float_local(sizes.fft_len);
    float *fft_i = alloc_aligned_float_local(sizes.fft_len);
    float *cobra_r = alloc_aligned_float_local(sizes.cobra_len);
    float *cobra_i = alloc_aligned_float_local(sizes.cobra_len);
    nufft1_plan *plan = nufft1_initialize_shared(m, n, 1, mode, tw_r, tw_i);
    nufft1_workspace *workspace = nufft1_create_workspace(plan, fft_r, fft_i, cobra_r, cobra_i);
    if (!plan || !workspace) return NAN;
    nufft1_precompute(workspace, x, m, 1.0);

    float *out_r = (float *)calloc((size_t)nout, sizeof(float));
    float *out_i = (float *)calloc((size_t)nout, sizeof(float));
    if (!out_r || !out_i) {
        free(out_r);
        free(out_i);
        nufft1_free_workspace(workspace);
        nufft1_free_plan(plan);
        free(tw_r);
        free(tw_i);
        free(fft_r);
        free(fft_i);
        free(cobra_r);
        free(cobra_i);
        return NAN;
    }

    nufft1_execute(workspace, yr, yi, out_r, out_i, 1);
    double total = 0.0;
    for (int r = 0; r < reps; ++r) {
        double t0 = now_seconds();
        nufft1_execute(workspace, yr, yi, out_r, out_i, 1);
        total += now_seconds() - t0;
    }

    free(out_r);
    free(out_i);
    nufft1_free_workspace(workspace);
    nufft1_free_plan(plan);
    free(tw_r);
    free(tw_i);
    free(fft_r);
    free(fft_i);
    free(cobra_r);
    free(cobra_i);
    return total / (double)reps;
}

int main(int argc, char **argv) {
    const char *path = argc > 1 ? argv[1] : "scaling.h";
    const int m = 512;
    const int n = 512;
    const int reps = 5;

    double *x = (double *)malloc((size_t)m * sizeof(double));
    float *yr = (float *)malloc((size_t)m * sizeof(float));
    float *yi = (float *)malloc((size_t)m * sizeof(float));
    if (!x || !yr || !yi) return 1;
    fill_inputs(x, yr, yi, m);

    double pre21 = time_precompute(NUFFT1_PSWF21, x, m, n, reps);
    double pre43 = time_precompute(NUFFT1_PSWF43, x, m, n, reps);
    double exe21 = time_execute(NUFFT1_PSWF21, x, yr, yi, m, n, n, reps);
    double exe43 = time_execute(NUFFT1_PSWF43, x, yr, yi, m, n, n + (n >> 1), reps);

    double beta21 = 2.5;
    double gamma21 = 10.0;
    double beta43 = 4.5;
    double gamma43 = 20.0;
    if (isfinite(exe21) && isfinite(exe43) && exe21 > 0.0) beta43 = fmax(1.0, beta21 * exe43 / exe21);
    if (isfinite(pre21) && isfinite(pre43) && pre21 > 0.0) gamma43 = fmax(1.0, gamma21 * pre43 / pre21);

    FILE *f = fopen(path, "w");
    if (!f) return 1;
    fprintf(f,
            "#define F_ALPHA %.10g\n"
            "#define F_PSWF21_BETA %.10g\n"
            "#define F_PSWF21_GAMMA %.10g\n"
            "#define F_PSWF43_BETA %.10g\n"
            "#define F_PSWF43_GAMMA %.10g\n",
            0.20, beta21, gamma21, beta43, gamma43);
    fclose(f);

    free(x);
    free(yr);
    free(yi);
    return 0;
}
