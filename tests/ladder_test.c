/* Regression tests for the two-level twiddle ladder (src/utils/ladder.h), the
 * vector trig kernels (src/utils/trig.h) and their use by the AoV sweeps.
 *
 * Build/run via tests/verify_ladder.sh.
 */
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils/aov.h"

static int failures;

static void check_close(const char *label, double got, double want, double tol) {
    double diff = fabs(got - want);
    if (!(diff <= tol)) {
        fprintf(stderr, "FAIL %s: got %.12g want %.12g (|diff|=%.3g > %.3g)\n", label, got, want, diff, tol);
        ++failures;
    }
}

static void check_true(const char *label, bool cond) {
    if (!cond) {
        fprintf(stderr, "FAIL %s\n", label);
        ++failures;
    }
}

/* Deterministic xorshift so the tests are reproducible without a library RNG. */
static uint64_t rng_state = UINT64_C(0x9e3779b97f4a7c15);

static double rnd_unit(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) / (double)(UINT64_C(1) << 53);
}

/* ------------------------------------------------------------------------- */
/* 1. Vector trig kernels must reproduce the scalar kernels exactly.          */
/* ------------------------------------------------------------------------- */
static void test_vector_trig(void) {
    double max_sin = 0.0;
    double max_cos = 0.0;

    for (int trial = 0; trial < 4000; ++trial) {
        VEC x;
        for (int lane = 0; lane < VEC_LEN; ++lane) {
            /* Mix small fractional phases with large ones, both signs. */
            double scale = (trial & 1) ? 1.0 : 1.0e5;
            x.values[lane] = (float)((rnd_unit() * 2.0 - 1.0) * scale);
        }
        VEC vs = sin2pif_vec(x);
        VEC vc = cos2pif_vec(x);
        for (int lane = 0; lane < VEC_LEN; ++lane) {
            double ds = fabs((double)vs.values[lane] - (double)sin2pif_tls(x.values[lane]));
            double dc = fabs((double)vc.values[lane] - (double)cos2pif_tls(x.values[lane]));
            if (ds > max_sin) max_sin = ds;
            if (dc > max_cos) max_cos = dc;
        }
    }
    printf("  vector trig vs scalar: max|dsin|=%.3e max|dcos|=%.3e\n", max_sin, max_cos);
    /* Identical arithmetic: allow only float-format noise. */
    check_close("vector sin parity", max_sin, 0.0, 1.0e-7);
    check_close("vector cos parity", max_cos, 0.0, 1.0e-7);
}

/* ------------------------------------------------------------------------- */
/* 2. Double phase fold and two-DVEC -> one-VEC concatenation.                */
/* ------------------------------------------------------------------------- */
static void test_phase_fold(void) {
    double max_fold = 0.0;
    double max_pair = 0.0;

    for (int trial = 0; trial < 4000; ++trial) {
        double phases[2 * DVEC_LEN];
        double scale = (trial & 1) ? 1.0 : 1.0e6;
        for (int lane = 0; lane < 2 * DVEC_LEN; ++lane) phases[lane] = (rnd_unit() * 2.0 - 1.0) * scale;

        DVEC lo, hi;
        for (int lane = 0; lane < DVEC_LEN; ++lane) {
            lo.values[lane] = phases[lane];
            hi.values[lane] = phases[DVEC_LEN + lane];
        }
        DVEC flo = aov_fold_phase_dvec(lo);
        DVEC fhi = aov_fold_phase_dvec(hi);
        for (int lane = 0; lane < DVEC_LEN; ++lane) {
            double want_lo = lo.values[lane] - (double)(long long)lo.values[lane];
            double want_hi = hi.values[lane] - (double)(long long)hi.values[lane];
            double d = fabs(flo.values[lane] - want_lo);
            if (d > max_fold) max_fold = d;
            d = fabs(fhi.values[lane] - want_hi);
            if (d > max_fold) max_fold = d;
        }

        VEC pair = aov_dvec_pair_to_vecf(flo, fhi);
        for (int lane = 0; lane < DVEC_LEN; ++lane) {
            double d = fabs((double)pair.values[lane] - flo.values[lane]);
            if (d > max_pair) max_pair = d;
            d = fabs((double)pair.values[DVEC_LEN + lane] - fhi.values[lane]);
            if (d > max_pair) max_pair = d;
        }
    }
    printf("  phase fold: max|state-exact|=%.3e  pair concat: max|diff|=%.3e\n", max_fold, max_pair);
    check_close("phase fold", max_fold, 0.0, 1.0e-12);
    check_close("pair concatenation", max_pair, 0.0, 1.0e-7);
}

/* ------------------------------------------------------------------------- */
/* Shared synthetic light curve.                                             */
/* ------------------------------------------------------------------------- */
enum { TEST_N = 500, TEST_PAD = 512 };

static double test_x[TEST_N];
static float test_y[TEST_N];
static float test_dy[TEST_N];
static float test_mult[TEST_N];

static void make_light_curve(double x0) {
    rng_state = UINT64_C(0x9e3779b97f4a7c15);
    for (int i = 0; i < TEST_N; ++i) {
        test_x[i] = x0 + 0.37 * (double)i + 3.0 * rnd_unit();
        test_y[i] = (float)(18.0 + 0.4 * sin(2.0 * M_PI * 0.13 * test_x[i]) + (rnd_unit() - 0.5));
        test_dy[i] = (float)(0.01 + 0.05 * rnd_unit());
        test_mult[i] = 1.0f / (test_dy[i] * test_dy[i]);
    }
}

/* ------------------------------------------------------------------------- */
/* 3. Ladder state must match the exact double-precision phase after any       */
/*    number of advances, with a bounded in-window drift.                     */
/* ------------------------------------------------------------------------- */
static void test_ladder_state(bool negative_times, bool use_mult) {
    const uint32_t block_len = 256;
    const size_t num_blocks = 4 * (size_t)MAX_TWIDDLE_REUSE;
    const double fmin = 0.0037;
    const double df = 1.0e-5;
    const uint32_t harmonics = 3;
    const float eps = 1.0e-6f;

    make_light_curve(negative_times ? -400.0 : 0.0);

    twiddle_ladder ladder;
    memset(&ladder, 0, sizeof(ladder));
    if (ladder_start(&ladder, test_x, TEST_N, TEST_PAD, harmonics, block_len, num_blocks, fmin, df, use_mult ? test_mult : NULL, test_dy, eps) != 0) {
        check_true("ladder_start", false);
        return;
    }

    /* The ladder advances one block at a time; check the state at every block so
     * anchors (block % MAX_TWIDDLE_REUSE == 0) and in-window blocks are both
     * covered, and that no third level is needed. */
    double max_anchor_err = 0.0;
    double max_drift = 0.0;

    for (size_t block = 0; block < num_blocks; ++block) {
        ladder_advance(&ladder, block);

        const double freq_block = fmin + df * (double)block_len * (double)block;
        for (uint32_t h = 0; h < harmonics; ++h) {
            const float *sr = NULL;
            const float *si = NULL;
            ladder_state(&ladder, h, &sr, &si);
            const double scale = (double)(h + 1U) * freq_block;
            double max_err = 0.0;
            for (int i = 0; i < TEST_N; ++i) {
                double phase = test_x[i] * scale;
                phase -= (double)(long long)phase;
                double amp = use_mult ? (double)test_mult[i] : 1.0 / ((double)test_dy[i] * (double)test_dy[i] + (double)eps);
                double wr = amp * cos(2.0 * M_PI * phase);
                double wi = amp * sin(2.0 * M_PI * phase);
                double er = fabs((double)sr[i] - wr);
                double ei = fabs((double)si[i] - wi);
                double err = fmax(er, ei) / fmax(1.0, fabs(amp));
                if (err > max_err) max_err = err;
            }
            if (block % (size_t)MAX_TWIDDLE_REUSE == 0U) {
                if (max_err > max_anchor_err) max_anchor_err = max_err;
            } else if (max_err > max_drift) {
                max_drift = max_err;
            }
        }
    }

    printf("  ladder state (%s times, %s): anchor max=%.3e in-window max=%.3e\n", negative_times ? "negative" : "positive",
           use_mult ? "mult" : "weight", max_anchor_err, max_drift);
    check_close("ladder anchor accuracy", max_anchor_err, 0.0, 2.0e-6);
    check_close("ladder in-window drift", max_drift, 0.0, 2.0e-5);

    ladder_free(&ladder);
}

/* ------------------------------------------------------------------------- */
/* 4. Merged block sums must reproduce brute-force double-precision sums, and */
/*    the solved AoV statistic must match the brute-force one.                */
/* ------------------------------------------------------------------------- */
static void *test_alloc(size_t bytes) {
    size_t rounded = (bytes + 63U) & ~(size_t)63U;
    if (rounded == 0U) rounded = 64U;
    return aligned_alloc(64, rounded);
}

static void test_block_sums(int degree) {
    const int grid_len = 512; /* NuFFT grid: block_len becomes nufft1_output_len() */
    const size_t num_blocks = 2;
    const double fmin = 0.0113;
    const double df = 4.0e-5;
    const int harmonics = 2 * degree;
    const float eps = 1.0e-6f;

    make_light_curve(0.0);

    /* Real NuFFT plan/workspace so the merged pass is validated end to end. */
    nufft1_external_sizes sizes = nufft1_external_sizes_for_plan(grid_len, NUFFT1_PSWF43);
    float *twiddle_r = test_alloc(sizes.twiddle_len * sizeof(float));
    float *twiddle_i = test_alloc(sizes.twiddle_len * sizeof(float));
    float *fft_r = test_alloc(sizes.fft_len * sizeof(float));
    float *fft_i = test_alloc(sizes.fft_len * sizeof(float));
    float *cobra_r = test_alloc(sizes.cobra_len * sizeof(float));
    float *cobra_i = test_alloc(sizes.cobra_len * sizeof(float));
    nufft1_plan *plan = NULL;
    nufft1_workspace *workspace = NULL;
    if (twiddle_r && twiddle_i) plan = nufft1_initialize_shared(TEST_N, grid_len, harmonics, NUFFT1_PSWF43, twiddle_r, twiddle_i);
    if (plan && fft_r && fft_i && cobra_r && cobra_i) workspace = nufft1_create_workspace(plan, fft_r, fft_i, cobra_r, cobra_i);
    check_true("nufft plan/workspace", plan != NULL && workspace != NULL);
    if (!workspace) return;

    nufft1_precompute(workspace, test_x, TEST_N, df);
    const uint32_t block_len = (uint32_t)nufft1_output_len(plan);
    const uint32_t nfreq = block_len + 137U; /* one full block + one partial */

    aov_reference_t ref;
    memset(&ref, 0, sizeof(ref));
    double ws = 0.0, wy = 0.0;
    for (int i = 0; i < TEST_N; ++i) {
        double w = 1.0 / ((double)test_dy[i] * (double)test_dy[i] + (double)eps);
        ws += w;
        wy += w * (double)test_y[i];
    }
    ref.ws = (float)ws;
    ref.ymean = (float)(wy / ws);
    double yws = 0.0, chi2_ref = 0.0;
    for (int i = 0; i < TEST_N; ++i) {
        double w = 1.0 / ((double)test_dy[i] * (double)test_dy[i] + (double)eps);
        double yc = (double)test_y[i] - (double)ref.ymean;
        yws += yc * w;
        chi2_ref += yc * yc * w;
    }
    ref.yws = (float)yws;
    ref.chi2_ref = (float)chi2_ref;

    const size_t trig_len = (size_t)(harmonics + 1) * block_len;
    const size_t ytrig_len = (size_t)(degree + 1) * block_len;
    float *ladder_Sw = calloc(trig_len, sizeof(float));
    float *ladder_Cw = calloc(trig_len, sizeof(float));
    float *ladder_Syw = calloc(ytrig_len, sizeof(float));
    float *ladder_Cyw = calloc(ytrig_len, sizeof(float));
    float *ref_Sw = calloc(trig_len, sizeof(float));
    float *ref_Cw = calloc(trig_len, sizeof(float));
    float *ref_Syw = calloc(ytrig_len, sizeof(float));
    float *ref_Cyw = calloc(ytrig_len, sizeof(float));
    float *pow_ladder = calloc(block_len, sizeof(float));
    float *cond_ladder = calloc(block_len, sizeof(float));
    float *pow_ref = calloc(block_len, sizeof(float));
    float *cond_ref = calloc(block_len, sizeof(float));
    float *tmp_r = calloc(block_len, sizeof(float));
    float *tmp_i = calloc(block_len, sizeof(float));
    bool allocs_ok = ladder_Sw && ladder_Cw && ladder_Syw && ladder_Cyw && ref_Sw && ref_Cw && ref_Syw && ref_Cyw && pow_ladder && cond_ladder &&
                     pow_ref && cond_ref && tmp_r && tmp_i;
    check_true("sums allocation", allocs_ok);
    if (!allocs_ok) return;

    twiddle_ladder ladder;
    memset(&ladder, 0, sizeof(ladder));
    check_true("ladder_start (block sums)",
               ladder_start(&ladder, test_x, TEST_N, TEST_PAD, (uint32_t)harmonics, block_len, num_blocks, fmin, df, NULL, test_dy, eps) == 0);

    buffer_t scratch;
    memset(&scratch, 0, sizeof(scratch));
    check_true("solver scratch", buffer_ensure_aov_solver_vec(&scratch, degree));

    double max_sum_err = 0.0;
    double max_r2_err = 0.0;

    for (size_t block = 0; block < num_blocks; ++block) {
        uint32_t base = (uint32_t)(block * (size_t)block_len);
        uint32_t count = block_len;
        if (base + count > nfreq) count = nfreq - base;
        const size_t stride = (size_t)count;

        if (block > 0U) ladder_advance(&ladder, block);
        ladder_block_sums(&ladder, workspace, test_y, ref.ymean, ladder_Sw, ladder_Cw, ladder_Syw, ladder_Cyw, tmp_r, tmp_i, count);

        /* Brute-force double-precision sums at the same frequencies. */
        memset(ref_Sw, 0, trig_len * sizeof(float));
        memset(ref_Cw, 0, trig_len * sizeof(float));
        memset(ref_Syw, 0, ytrig_len * sizeof(float));
        memset(ref_Cyw, 0, ytrig_len * sizeof(float));
        for (uint32_t i = 0; i < count; ++i) {
            const double freq = fmin + (double)(base + i) * df;
            for (int j = 0; j < TEST_N; ++j) {
                double w = 1.0 / ((double)test_dy[j] * (double)test_dy[j] + (double)eps);
                double yc = (double)test_y[j] - (double)ref.ymean;
                for (int q = 1; q <= harmonics; ++q) {
                    double phase = test_x[j] * (double)q * freq;
                    phase -= (double)(long long)phase;
                    double ang = 2.0 * M_PI * phase;
                    double c = cos(ang), s = sin(ang);
                    ref_Cw[(size_t)q * stride + i] += (float)(w * c);
                    ref_Sw[(size_t)q * stride + i] += (float)(w * s);
                    if (q <= degree) {
                        ref_Cyw[(size_t)q * stride + i] += (float)(w * yc * c);
                        ref_Syw[(size_t)q * stride + i] += (float)(w * yc * s);
                    }
                }
            }
            for (int q = 1; q <= harmonics; ++q) {
                double dl = fabs(hypot((double)ladder_Cw[(size_t)q * stride + i], (double)ladder_Sw[(size_t)q * stride + i]) -
                                 hypot((double)ref_Cw[(size_t)q * stride + i], (double)ref_Sw[(size_t)q * stride + i])) /
                            (double)ws;
                if (dl > max_sum_err) max_sum_err = dl;
                if (q <= degree) {
                    double dy = fabs(hypot((double)ladder_Cyw[(size_t)q * stride + i], (double)ladder_Syw[(size_t)q * stride + i]) -
                                     hypot((double)ref_Cyw[(size_t)q * stride + i], (double)ref_Syw[(size_t)q * stride + i])) /
                                (double)ws;
                    if (dy > max_sum_err) max_sum_err = dy;
                }
            }
        }
        /* The q = 0 rows are constants over the whole run. */
        for (uint32_t i = 0; i < count; ++i) {
            ref_Cw[i] = ref.ws;
            ref_Sw[i] = 0.0f;
            ref_Cyw[i] = ref.yws;
            ref_Syw[i] = 0.0f;
            ladder_Cw[i] = ref.ws;
            ladder_Sw[i] = 0.0f;
            ladder_Cyw[i] = ref.yws;
            ladder_Syw[i] = 0.0f;
        }

        if (degree == 1) {
            aov_gls_impl(ladder_Sw, ladder_Cw, ladder_Syw, ladder_Cyw, (int)count, ref.ws, ref.yws, ref.chi2_ref, pow_ladder, cond_ladder);
            aov_gls_impl(ref_Sw, ref_Cw, ref_Syw, ref_Cyw, (int)count, ref.ws, ref.yws, ref.chi2_ref, pow_ref, cond_ref);
        } else {
            if (aov_solve_periodogram_vec(ladder_Sw, ladder_Cw, ladder_Syw, ladder_Cyw, (int)count, degree, ref.chi2_ref, pow_ladder, cond_ladder,
                                          &scratch) != 0 ||
                aov_solve_periodogram_vec(ref_Sw, ref_Cw, ref_Syw, ref_Cyw, (int)count, degree, ref.chi2_ref, pow_ref, cond_ref, &scratch) != 0) {
                check_true("periodogram solve", false);
                break;
            }
        }

        for (uint32_t i = 0; i < count; ++i) {
            if (!float_is_finite_bits(pow_ladder[i]) || !float_is_finite_bits(pow_ref[i])) continue;
            double d = fabs((double)pow_ladder[i] - (double)pow_ref[i]);
            if (d > max_r2_err) max_r2_err = d;
        }
    }

    printf("  block sums d=%d: max relative sum error=%.3e  max |dr2|=%.3e\n", degree, max_sum_err, max_r2_err);
    /* The NuFFT itself is accurate to ~2e-5 relative, so that bounds the sums. */
    check_close("block sums vs direct", max_sum_err, 0.0, 1.0e-4);
    check_close("solved r2 vs direct", max_r2_err, 0.0, 1.0e-5);

    free(tmp_r);
    free(tmp_i);
    free(ladder_Sw);
    free(ladder_Cw);
    free(ladder_Syw);
    free(ladder_Cyw);
    free(ref_Sw);
    free(ref_Cw);
    free(ref_Syw);
    free(ref_Cyw);
    free(pow_ladder);
    free(cond_ladder);
    free(pow_ref);
    free(cond_ref);
    free(scratch.aovSolverVec);
    ladder_free(&ladder);
    nufft1_free_workspace(workspace);
    nufft1_free_plan(plan);
    free(twiddle_r);
    free(twiddle_i);
    free(fft_r);
    free(fft_i);
    free(cobra_r);
    free(cobra_i);
}

int main(void) {
    printf("ladder_test: two-level twiddle ladder and vector trig\n");

    printf("[1] vector trig parity\n");
    test_vector_trig();

    printf("[2] double phase fold\n");
    test_phase_fold();

    printf("[3] ladder state accuracy\n");
    test_ladder_state(false, false);
    test_ladder_state(true, false);
    test_ladder_state(false, true);

    printf("[4] merged block sums vs direct sums\n");
    test_block_sums(1);
    test_block_sums(3);

    if (failures != 0) {
        fprintf(stderr, "ladder_test: %d FAILURE(S)\n", failures);
        return 1;
    }
    printf("ladder_test: PASS\n");
    return 0;
}
