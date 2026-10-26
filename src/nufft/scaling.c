#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../utils/compat.h"
#include "nufft1.h"

/* Build the calibration helper as a single translation unit so it can call the
 * static COBRA/COBRAVO kernels directly.  nufft1.c must not include the
 * generated scaling.h while this TU is producing it. */
#ifndef NUFFT1_SCALING_BOOTSTRAP
#    define NUFFT1_SCALING_BOOTSTRAP 1
#endif
#include "nufft1.c"

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
    double best = 1e30;
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
        double dt = now_seconds() - t0;
        if (dt < best) best = dt;
        nufft1_free_workspace(workspace);
        nufft1_free_plan(plan);
        free(tw_r);
        free(tw_i);
        free(fft_r);
        free(fft_i);
        free(cobra_r);
        free(cobra_i);
    }
    return best;
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
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        double t0 = now_seconds();
        nufft1_execute(workspace, yr, yi, out_r, out_i, 1);
        double dt = now_seconds() - t0;
        if (dt < best) best = dt;
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
    return best;
}

static void fill_fft_inputs(float *real, float *imag, size_t n, unsigned int seed) {
    srand(seed);
    for (size_t k = 0; k < n; ++k) {
        real[k] = (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
        imag[k] = (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
    }
}

/* Time one FFT kernel on a cold-but-resident input.  The copy is outside the
 * timed region so repeated transforms do not overflow the float range. */
static double fft_time(int use_fused, uint32_t n, int reps) {
    float *real = alloc_aligned_float_local(n);
    float *imag = alloc_aligned_float_local(n);
    float *orig_real = alloc_aligned_float_local(n);
    float *orig_imag = alloc_aligned_float_local(n);
    float *tw_r = alloc_aligned_float_local(n);
    float *tw_i = alloc_aligned_float_local(n);
    if (!real || !imag || !orig_real || !orig_imag || !tw_r || !tw_i) {
        free(real);
        free(imag);
        free(orig_real);
        free(orig_imag);
        free(tw_r);
        free(tw_i);
        return NAN;
    }
    fill_fft_inputs(orig_real, orig_imag, n, 1234U);
    if (use_fused)
        generate_fft_buffer_fused_r2(n, tw_r, tw_i);
    else
        generate_fft_buffer(n, tw_r, tw_i);

    double best = 1e30;
    for (int rep = 0; rep <= reps; ++rep) {
        memcpy(real, orig_real, (size_t)n * sizeof(float));
        memcpy(imag, orig_imag, (size_t)n * sizeof(float));
        double t0 = now_seconds();
        if (use_fused)
            sande_tukey_in_place_fused_r2(real, imag, tw_r, tw_i, n);
        else
            sande_tukey_in_place(real, imag, tw_r, tw_i, n);
        double dt = now_seconds() - t0;
        if (rep > 0 && dt < best) best = dt;
    }
    free(real);
    free(imag);
    free(orig_real);
    free(orig_imag);
    free(tw_r);
    free(tw_i);
    return best;
}

static void naive_dft(const float *real, const float *imag, float *out_r, float *out_i, size_t n) {
    for (size_t k = 0; k < n; ++k) {
        double sr = 0.0, si = 0.0;
        for (size_t t = 0; t < n; ++t) {
            double a = 2.0 * M_PI * (double)((k * t) % n) / (double)n;
            double c = cos(a), s = sin(a);
            sr += (double)real[t] * c - (double)imag[t] * s;
            si += (double)real[t] * s + (double)imag[t] * c;
        }
        out_r[k] = (float)sr;
        out_i[k] = (float)si;
    }
}

static int fft_self_test(void) {
    static const uint32_t orders[] = {10U, 13U, 14U, 15U, 16U, 17U};
    for (size_t t = 0; t < sizeof(orders) / sizeof(orders[0]); ++t) {
        uint32_t order = orders[t];
        uint32_t n = 1U << order;
        float *ref_r = alloc_aligned_float_local(n);
        float *ref_i = alloc_aligned_float_local(n);
        float *got_r = alloc_aligned_float_local(n);
        float *got_i = alloc_aligned_float_local(n);
        float *tw_old_r = alloc_aligned_float_local(n);
        float *tw_old_i = alloc_aligned_float_local(n);
        float *tw_fused_r = alloc_aligned_float_local(n);
        float *tw_fused_i = alloc_aligned_float_local(n);
        if (!ref_r || !ref_i || !got_r || !got_i || !tw_old_r || !tw_old_i || !tw_fused_r || !tw_fused_i) {
            free(ref_r);
            free(ref_i);
            free(got_r);
            free(got_i);
            free(tw_old_r);
            free(tw_old_i);
            free(tw_fused_r);
            free(tw_fused_i);
            return 1;
        }
        fill_fft_inputs(ref_r, ref_i, n, 900U + order);
        memcpy(got_r, ref_r, (size_t)n * sizeof(float));
        memcpy(got_i, ref_i, (size_t)n * sizeof(float));
        generate_fft_buffer(n, tw_old_r, tw_old_i);
        generate_fft_buffer_fused_r2(n, tw_fused_r, tw_fused_i);
        sande_tukey_in_place(ref_r, ref_i, tw_old_r, tw_old_i, n);
        sande_tukey_in_place_fused_r2(got_r, got_i, tw_fused_r, tw_fused_i, n);
        table_shuffle(ref_r, ref_i, order);
        table_shuffle(got_r, got_i, order);

        double max_abs = 0.0, max_rel = 0.0;
        for (uint32_t k = 0; k < n; ++k) {
            double dr = fabs((double)got_r[k] - (double)ref_r[k]);
            double di = fabs((double)got_i[k] - (double)ref_i[k]);
            double mag = sqrt((double)ref_r[k] * (double)ref_r[k] + (double)ref_i[k] * (double)ref_i[k]);
            if (dr > max_abs) max_abs = dr;
            if (di > max_abs) max_abs = di;
            if (mag > 1e-9) {
                double rel = sqrt(dr * dr + di * di) / mag;
                if (rel > max_rel) max_rel = rel;
            }
        }
        if (!(max_abs < 5e-3 && max_rel < 5e-4)) {
            fprintf(stderr, "scaling: fused FFT mismatch at N=%u max_abs=%.3e max_rel=%.3e\n", n, max_abs, max_rel);
            free(ref_r);
            free(ref_i);
            free(got_r);
            free(got_i);
            free(tw_old_r);
            free(tw_old_i);
            free(tw_fused_r);
            free(tw_fused_i);
            return 1;
        }
        if (order == 10U) {
            /* Reference the fused result against a direct DFT as well as the
             * legacy kernel, so a shared bug cannot pass the self-test. */
            float *in_r = alloc_aligned_float_local(n);
            float *in_i = alloc_aligned_float_local(n);
            float *out_r = alloc_aligned_float_local(n);
            float *out_i = alloc_aligned_float_local(n);
            if (!in_r || !in_i || !out_r || !out_i) {
                free(in_r);
                free(in_i);
                free(out_r);
                free(out_i);
                free(ref_r);
                free(ref_i);
                free(got_r);
                free(got_i);
                free(tw_old_r);
                free(tw_old_i);
                free(tw_fused_r);
                free(tw_fused_i);
                return 1;
            }
            fill_fft_inputs(in_r, in_i, n, 777U);
            naive_dft(in_r, in_i, out_r, out_i, n);
            /* Re-run the fused transform on the same input; the previous
             * got_* buffers were consumed by table_shuffle. */
            memcpy(got_r, in_r, (size_t)n * sizeof(float));
            memcpy(got_i, in_i, (size_t)n * sizeof(float));
            sande_tukey_in_place_fused_r2(got_r, got_i, tw_fused_r, tw_fused_i, n);
            table_shuffle(got_r, got_i, order);
            double dft_abs = 0.0, dft_rel = 0.0;
            for (uint32_t k = 0; k < n; ++k) {
                double dr = fabs((double)got_r[k] - (double)out_r[k]);
                double di = fabs((double)got_i[k] - (double)out_i[k]);
                double mag = sqrt((double)out_r[k] * (double)out_r[k] + (double)out_i[k] * (double)out_i[k]);
                if (dr > dft_abs) dft_abs = dr;
                if (di > dft_abs) dft_abs = di;
                if (mag > 1e-9) {
                    double rel = sqrt(dr * dr + di * di) / mag;
                    if (rel > dft_rel) dft_rel = rel;
                }
            }
            if (!(dft_abs < 2e-3 && dft_rel < 2e-4)) {
                fprintf(stderr, "scaling: fused FFT vs DFT mismatch at N=%u max_abs=%.3e max_rel=%.3e\n", n, dft_abs, dft_rel);
                free(in_r);
                free(in_i);
                free(out_r);
                free(out_i);
                free(ref_r);
                free(ref_i);
                free(got_r);
                free(got_i);
                free(tw_old_r);
                free(tw_old_i);
                free(tw_fused_r);
                free(tw_fused_i);
                return 1;
            }
            free(in_r);
            free(in_i);
            free(out_r);
            free(out_i);
        }
        printf("scaling: fused FFT N=%u max_rel=%.3e OK\n", n, max_rel);
        free(ref_r);
        free(ref_i);
        free(got_r);
        free(got_i);
        free(tw_old_r);
        free(tw_old_i);
        free(tw_fused_r);
        free(tw_fused_i);
    }
    return 0;
}

static int fit_power_alpha(double *alpha_out, double *r2_out) {
    /* Skip the small sizes: their fixed overhead biases the log-log slope
     * and made the build-time alpha visibly noisy. */
    static const uint32_t ns[] = {16384U, 32768U, 65536U, 131072U, 262144U};
    const size_t count = sizeof(ns) / sizeof(ns[0]);
    double xs[sizeof(ns) / sizeof(ns[0])];
    double ys[sizeof(ns) / sizeof(ns[0])];
    for (size_t i = 0; i < count; ++i) {
        double t = fft_time(1, ns[i], 20);
        if (!isfinite(t) || t <= 0.0) return 1;
        xs[i] = log((double)ns[i]);
        ys[i] = log(t);
    }
    double mean_x = 0.0, mean_y = 0.0;
    for (size_t i = 0; i < count; ++i) {
        mean_x += xs[i];
        mean_y += ys[i];
    }
    mean_x /= (double)count;
    mean_y /= (double)count;
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < count; ++i) {
        num += (xs[i] - mean_x) * (ys[i] - mean_y);
        den += (xs[i] - mean_x) * (xs[i] - mean_x);
    }
    if (den <= 0.0) return 1;
    double slope = num / den;
    double intercept = mean_y - slope * mean_x;
    double ss_tot = 0.0, ss_res = 0.0;
    for (size_t i = 0; i < count; ++i) {
        double pred = intercept + slope * xs[i];
        ss_tot += (ys[i] - mean_y) * (ys[i] - mean_y);
        ss_res += (ys[i] - pred) * (ys[i] - pred);
    }
    double r2 = ss_tot > 0.0 ? 1.0 - ss_res / ss_tot : 0.0;
    double alpha = slope - 1.0;
    if (!isfinite(alpha) || !isfinite(r2) || r2 < 0.90) return 1;
    if (alpha < 0.10) alpha = 0.10;
    if (alpha > 0.30) alpha = 0.30;
    *alpha_out = alpha;
    *r2_out = r2;
    return 0;
}

static int bitrev_self_test(void) {
    static const uint32_t orders[] = {13U, 14U, 16U};
    for (size_t t = 0; t < sizeof(orders) / sizeof(orders[0]); ++t) {
        uint32_t order = orders[t];
        size_t n = (size_t)1U << order;
        float *ref_r = alloc_aligned_float_local(n);
        float *ref_i = alloc_aligned_float_local(n);
        float *got_r = alloc_aligned_float_local(n);
        float *got_i = alloc_aligned_float_local(n);
        float *buf_a_r = alloc_aligned_float_local(2U * (size_t)BLOCK_ELEMS);
        float *buf_a_i = alloc_aligned_float_local(2U * (size_t)BLOCK_ELEMS);
        float *buf_b_r = alloc_aligned_float_local((size_t)BLOCK_ELEMS);
        float *buf_b_i = alloc_aligned_float_local((size_t)BLOCK_ELEMS);
        if (!ref_r || !ref_i || !got_r || !got_i || !buf_a_r || !buf_a_i || !buf_b_r || !buf_b_i) {
            free(ref_r);
            free(ref_i);
            free(got_r);
            free(got_i);
            free(buf_a_r);
            free(buf_a_i);
            free(buf_b_r);
            free(buf_b_i);
            return 1;
        }
        for (size_t k = 0; k < n; ++k) {
            ref_r[k] = got_r[k] = (float)k + 0.25f;
            ref_i[k] = got_i[k] = -(float)k - 0.5f;
        }
        table_shuffle(ref_r, ref_i, order);

        cobra_shuffle(got_r, got_i, order, buf_a_r, buf_a_i);
        if (memcmp(ref_r, got_r, n * sizeof(float)) != 0 || memcmp(ref_i, got_i, n * sizeof(float)) != 0) {
            fprintf(stderr, "scaling: COBRA bit-reversal mismatch at N=%zu\n", n);
            return 1;
        }

        for (size_t k = 0; k < n; ++k) {
            got_r[k] = (float)k + 0.25f;
            got_i[k] = -(float)k - 0.5f;
        }
        cobravo_init();
        cobravo_shuffle(got_r, got_i, order, buf_a_r, buf_a_i, buf_b_r, buf_b_i);
        if (memcmp(ref_r, got_r, n * sizeof(float)) != 0 || memcmp(ref_i, got_i, n * sizeof(float)) != 0) {
            fprintf(stderr, "scaling: COBRAVO bit-reversal mismatch at N=%zu\n", n);
            return 1;
        }
        free(ref_r);
        free(ref_i);
        free(got_r);
        free(got_i);
        free(buf_a_r);
        free(buf_a_i);
        free(buf_b_r);
        free(buf_b_i);
    }
    return 0;
}

static double bitrev_time(int use_cobravo, uint32_t order, int reps) {
    size_t n = (size_t)1U << order;
    float *real = alloc_aligned_float_local(n);
    float *imag = alloc_aligned_float_local(n);
    float *buf_a_r = alloc_aligned_float_local(2U * (size_t)BLOCK_ELEMS);
    float *buf_a_i = alloc_aligned_float_local(2U * (size_t)BLOCK_ELEMS);
    float *buf_b_r = alloc_aligned_float_local((size_t)BLOCK_ELEMS);
    float *buf_b_i = alloc_aligned_float_local((size_t)BLOCK_ELEMS);
    if (!real || !imag || !buf_a_r || !buf_a_i || !buf_b_r || !buf_b_i) {
        free(real);
        free(imag);
        free(buf_a_r);
        free(buf_a_i);
        free(buf_b_r);
        free(buf_b_i);
        return NAN;
    }
    for (size_t k = 0; k < n; ++k) {
        real[k] = (float)k;
        imag[k] = -(float)k;
    }
    cobravo_init();
    for (int w = 0; w < 3; ++w) {
        if (use_cobravo)
            cobravo_shuffle(real, imag, order, buf_a_r, buf_a_i, buf_b_r, buf_b_i);
        else
            cobra_shuffle(real, imag, order, buf_a_r, buf_a_i);
    }
    double total = 0.0;
    for (int r = 0; r < reps; ++r) {
        double t0 = now_seconds();
        if (use_cobravo)
            cobravo_shuffle(real, imag, order, buf_a_r, buf_a_i, buf_b_r, buf_b_i);
        else
            cobra_shuffle(real, imag, order, buf_a_r, buf_a_i);
        total += now_seconds() - t0;
    }
    free(real);
    free(imag);
    free(buf_a_r);
    free(buf_a_i);
    free(buf_b_r);
    free(buf_b_i);
    return total / (double)reps;
}

static int bitrev_benchmark(int *use_cobravo) {
    const uint32_t order = 16U;
    const int reps = 25;
    double cobra = bitrev_time(0, order, reps);
    double cobravo = bitrev_time(1, order, reps);
    if (!isfinite(cobra) || !isfinite(cobravo)) return 1;
    /* Run each kernel a second time to reduce sensitivity to scheduler noise. */
    double cobra2 = bitrev_time(0, order, reps);
    double cobravo2 = bitrev_time(1, order, reps);
    if (isfinite(cobra2)) cobra = fmin(cobra, cobra2);
    if (isfinite(cobravo2)) cobravo = fmin(cobravo, cobravo2);
    *use_cobravo = (cobravo <= cobra) ? 1 : 0;
    printf("scaling: bit-reversal N=%u COBRA %.3f us, COBRAVO %.3f us, %.2fx -> %s\n", 1U << order, cobra * 1e6, cobravo * 1e6,
           cobra > 0.0 ? cobra / cobravo : 0.0, *use_cobravo ? "COBRAVO" : "COBRA");
    fflush(stdout);
    return 0;
}

int main(int argc, char **argv) {
    const char *path = argc > 1 ? argv[1] : "scaling.h";
    const int m = 512;
    const int n = 512;
    const int reps = 30;

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
    printf("scaling: pre21=%.3f us pre43=%.3f us exe21=%.3f us exe43=%.3f us beta43=%.6f gamma43=%.6f\n", pre21 * 1e6, pre43 * 1e6, exe21 * 1e6, exe43 * 1e6,
           beta43, gamma43);
    fflush(stdout);

    if (bitrev_self_test() != 0) return 1;
    if (fft_self_test() != 0) return 1;
    int use_cobravo = 1;
    if (bitrev_benchmark(&use_cobravo) != 0) return 1;

    double alpha = 0.20;
    double alpha_r2 = 0.0;
    if (fit_power_alpha(&alpha, &alpha_r2) == 0) {
        printf("scaling: fitted power alpha=%.6f (R2=%.6f)\n", alpha, alpha_r2);
    } else {
        printf("scaling: power-alpha fit failed; using fallback %.6f\n", alpha);
    }
    fflush(stdout);

    double fft_legacy = fft_time(0, 65536U, 8);
    double fft_fused = fft_time(1, 65536U, 8);
    if (isfinite(fft_legacy) && isfinite(fft_fused) && fft_fused > 0.0) {
        printf("scaling: FFT N=65536 legacy %.3f us, fused %.3f us, %.2fx\n", fft_legacy * 1e6, fft_fused * 1e6, fft_legacy / fft_fused);
        fflush(stdout);
    }

    FILE *f = fopen(path, "w");
    if (!f) return 1;
    fprintf(f,
            "#define F_ALPHA %.10g\n"
            "#define F_PSWF21_BETA %.10g\n"
            "#define F_PSWF21_GAMMA %.10g\n"
            "#define F_PSWF43_BETA %.10g\n"
            "#define F_PSWF43_GAMMA %.10g\n"
            "#define F_BITREV_COBRAVO %d\n",
            alpha, beta21, gamma21, beta43, gamma43, use_cobravo);
    fclose(f);

    free(x);
    free(yr);
    free(yi);
    return 0;
}
