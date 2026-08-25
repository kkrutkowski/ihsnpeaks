#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "../src/utils/common.h"
#include "../src/utils/readout.h"

// Solve A * theta = b via Gaussian elimination with partial pivoting (double precision)
static bool solve_linear_system(int m, double* A, double* b, double* theta) {
    // Augmented matrix [A | b]
    double* M = (double*)malloc((size_t)m * (size_t)(m + 1) * sizeof(double));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            M[i * (m + 1) + j] = A[i * m + j];
        }
        M[i * (m + 1) + m] = b[i];
    }

    for (int i = 0; i < m; i++) {
        // Pivot
        int max_row = i;
        double max_val = fabs(M[i * (m + 1) + i]);
        for (int r = i + 1; r < m; r++) {
            if (fabs(M[r * (m + 1) + i]) > max_val) {
                max_val = fabs(M[r * (m + 1) + i]);
                max_row = r;
            }
        }
        if (max_val < 1e-15) {
            free(M);
            return false;
        }
        if (max_row != i) {
            for (int c = i; c <= m; c++) {
                double tmp = M[i * (m + 1) + c];
                M[i * (m + 1) + c] = M[max_row * (m + 1) + c];
                M[max_row * (m + 1) + c] = tmp;
            }
        }

        double pivot = M[i * (m + 1) + i];
        for (int r = i + 1; r < m; r++) {
            double factor = M[r * (m + 1) + i] / pivot;
            for (int c = i; c <= m; c++) {
                M[r * (m + 1) + c] -= factor * M[i * (m + 1) + c];
            }
        }
    }

    // Back substitution
    for (int i = m - 1; i >= 0; i--) {
        double sum = M[i * (m + 1) + m];
        for (int j = i + 1; j < m; j++) {
            sum -= M[i * (m + 1) + j] * theta[j];
        }
        theta[i] = sum / M[i * (m + 1) + i];
    }

    free(M);
    return true;
}

// Direct weighted least squares fit of the exact trig polynomial + 1/2 harmonic basis
static void direct_trig_fit(const double* x, const float* y, const float* dy, uint32_t n, int degree, double epsilon, float* y_detrended_out, double* magnitude_out) {
    double span = x[n - 1] - x[0];
    double ws = 0.0;
    double wysum = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        double w = 1.0;
        if (dy) {
            double dy_val = (double)dy[i];
            if (dy_val > 0.0 && isfinite(dy_val)) {
                w = 1.0 / (dy_val * dy_val + epsilon);
            } else {
                w = 1.0 / (999.9 * 999.9 + epsilon);
            }
        }
        ws += w;
        wysum += w * (double)y[i];
    }
    double y_mean = ws > 0.0 ? wysum / ws : 0.0;
    *magnitude_out = y_mean;

    int max_supported_degree = ((int)n - 2) / 4;
    if (max_supported_degree < 0) max_supported_degree = 0;
    int eff_degree = degree < max_supported_degree ? degree : max_supported_degree;

    if (span <= 0.0 || n < 6 || eff_degree == 0) {
        for (uint32_t i = 0; i < n; i++) {
            y_detrended_out[i] = (float)((double)y[i] - y_mean);
        }
        return;
    }

    // Number of basis columns: 1 (constant) + 2 * eff_degree (harmonics 1..eff_degree at f_base = 1/(2*span))
    int nbasis = 1 + 2 * eff_degree;
    double* A_basis = (double*)malloc((size_t)n * (size_t)nbasis * sizeof(double));
    double f_base = 0.5 / span;
    double t_mid = (x[0] + x[n - 1]) / 2.0;

    for (uint32_t i = 0; i < n; i++) {
        double phase = (x[i] - t_mid) * f_base;
        A_basis[i * nbasis + 0] = 1.0; // constant
        for (int h = 1; h <= eff_degree; h++) {
            A_basis[i * nbasis + 1 + 2 * (h - 1)] = sin(2.0 * M_PI * (double)h * phase);
            A_basis[i * nbasis + 2 + 2 * (h - 1)] = cos(2.0 * M_PI * (double)h * phase);
        }
    }

    // Normal equations: (A^T W A) theta = A^T W y
    double* ATA = (double*)calloc((size_t)nbasis * (size_t)nbasis, sizeof(double));
    double* ATy = (double*)calloc((size_t)nbasis, sizeof(double));

    for (uint32_t i = 0; i < n; i++) {
        double w = 1.0;
        if (dy) {
            double dy_val = (double)dy[i];
            if (dy_val > 0.0 && isfinite(dy_val)) {
                w = 1.0 / (dy_val * dy_val + epsilon);
            } else {
                w = 1.0 / (999.9 * 999.9 + epsilon);
            }
        }
        for (int j = 0; j < nbasis; j++) {
            double a_ij = A_basis[i * nbasis + j];
            ATy[j] += w * a_ij * (double)y[i];
            for (int k = 0; k < nbasis; k++) {
                ATA[j * nbasis + k] += w * a_ij * A_basis[i * nbasis + k];
            }
        }
    }

    double* theta = (double*)calloc((size_t)nbasis, sizeof(double));
    bool solved = solve_linear_system(nbasis, ATA, ATy, theta);
    assert(solved);

    for (uint32_t i = 0; i < n; i++) {
        double fit_val = 0.0;
        for (int j = 0; j < nbasis; j++) {
            fit_val += theta[j] * A_basis[i * nbasis + j];
        }
        y_detrended_out[i] = (float)((double)y[i] - fit_val);
    }

    free(A_basis);
    free(ATA);
    free(ATy);
    free(theta);
}

static void run_test_case(const char* test_name, uint32_t n, int degree, bool non_uniform_grid, bool weighted) {
    buffer_t buf = {0};
    buf.n = n;
    buf.x = (double*)malloc(n * sizeof(double));
    buf.y = (float*)malloc(n * sizeof(float));
    buf.dy = (float*)malloc(n * sizeof(float));

    double t_span = 100.0;
    for (uint32_t i = 0; i < n; i++) {
        if (non_uniform_grid) {
            double frac = (double)i / (double)(n - 1);
            buf.x[i] = 10.0 + t_span * (frac * frac + 0.1 * sin(5.0 * frac));
        } else {
            buf.x[i] = 10.0 + t_span * ((double)i / (double)(n - 1));
        }
        if (weighted) {
            buf.dy[i] = 0.01f + 0.05f * (float)((i % 7) + 1) / 7.0f;
        } else {
            buf.dy[i] = 0.05f;
        }

        // True signal = baseline + trend + oscillation + noise
        double t = buf.x[i] - buf.x[0];
        double trend = 14.5 + 0.01 * t + 0.5 * sin(M_PI * (t - t_span/2.0)/t_span) + 0.3 * cos(2.0 * M_PI * (t - t_span/2.0)/t_span);
        double signal = 0.1 * sin(2.0 * M_PI * 1.5 * t); // period = 1/1.5
        buf.y[i] = (float)(trend + signal + 0.002 * sin(123.4 * (double)i));
    }

    // Save copy for direct fit
    double* orig_x = (double*)malloc(n * sizeof(double));
    float* orig_y = (float*)malloc(n * sizeof(float));
    float* orig_dy = (float*)malloc(n * sizeof(float));
    memcpy(orig_x, buf.x, n * sizeof(double));
    memcpy(orig_y, buf.y, n * sizeof(float));
    memcpy(orig_dy, buf.dy, n * sizeof(float));

    double epsilon = 0.001;
    detrend_buffer_szego(&buf, epsilon, degree);

    float* ref_y_detrended = (float*)malloc(n * sizeof(float));
    double ref_magnitude = 0.0;
    direct_trig_fit(orig_x, orig_y, orig_dy, n, degree, epsilon, ref_y_detrended, &ref_magnitude);

    // Check magnitude
    double mag_err = fabs((double)buf.magnitude - ref_magnitude);
    assert(mag_err < 1e-4);

    // Check detrended values agreement with exact least squares
    double max_err = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        double err = fabs((double)buf.y[i] - (double)ref_y_detrended[i]);
        if (err > max_err) max_err = err;
    }

    printf("[PASS] %s (n=%u, degree=%d, max_diff=%.3e)\n", test_name, n, degree, max_err);
    assert(max_err < 1e-4); // Float precision match

    free(buf.x); free(buf.y); free(buf.dy);
    free(orig_x); free(orig_y); free(orig_dy);
    free(ref_y_detrended);
}

int main(void) {
    printf("Running Detrending Test Suite...\n\n");

    // 1. Basic degree tests (uniform, unweighted)
    run_test_case("Uniform degree=0 (constant)", 100, 0, false, false);
    run_test_case("Uniform degree=1 (0.5 f0)", 100, 1, false, false);
    run_test_case("Uniform degree=2 (0.5 f0 + 1.0 f0)", 100, 2, false, false);
    run_test_case("Uniform degree=3 (0.5 f0 + 1.0 f0 + 1.5 f0, default)", 100, 3, false, false);
    run_test_case("Uniform degree=5 (0.5 .. 2.5 f0)", 200, 5, false, false);

    // 2. Non-uniform sampling & weighted tests
    run_test_case("Non-uniform weighted degree=0", 150, 0, true, true);
    run_test_case("Non-uniform weighted degree=1", 150, 1, true, true);
    run_test_case("Non-uniform weighted degree=2", 150, 2, true, true);
    run_test_case("Non-uniform weighted degree=3", 150, 3, true, true);
    run_test_case("Non-uniform weighted degree=4", 150, 4, true, true);

    // 3. Edge case clamping tests
    run_test_case("Edge case n=5 requested degree=3 (clamps to d=0)", 5, 3, false, false);
    run_test_case("Edge case n=7 requested degree=3 (clamps to d=1)", 7, 3, false, false);
    run_test_case("Edge case n=11 requested degree=3 (clamps to d=2)", 11, 3, false, false);
    run_test_case("Edge case n=14 requested degree=3 (supports d=3)", 14, 3, false, false);

    printf("\nAll 14 tests PASSED successfully!\n");
    return 0;
}
