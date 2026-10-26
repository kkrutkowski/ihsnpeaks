/* Regression tests for the d=1 condition/Cybenko path, the forward-only
 * Levinson solver, and the vectorized IHS Wilkie correction.
 *
 * Build/run via tests/verify_core.sh.
 */
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils/readout.h"

/* process.h-only helpers referenced by aov.h's streaming code.  The tests
 * below only exercise aov_gls_impl() and the Levinson solver, so these stubs
 * are never called; they only satisfy the header. */
static inline void fill_twiddle_ladder(buffer_t *buffer, double fstep,
                                       double freq_factor) {
  (void)buffer;
  (void)fstep;
  (void)freq_factor;
}
static inline void compute_work_at_block(buffer_t *buffer, size_t block_idx) {
  (void)buffer;
  (void)block_idx;
}
static inline void advance_work_ladder(buffer_t *buffer, size_t next_block) {
  (void)buffer;
  (void)next_block;
}
static inline void appendFreq(double freq, float magnitude, int precision,
                              bool outputPeriod, sds *buffer,
                              char *stringBuff) {
  (void)freq;
  (void)magnitude;
  (void)precision;
  (void)outputPeriod;
  (void)buffer;
  (void)stringBuff;
}

#include "utils/aov.h"

static int failures;

static void check_close(const char *label, double got, double want,
                        double tol) {
  double diff = fabs(got - want);
  if (!(diff <= tol)) {
    fprintf(stderr, "FAIL %s: got %.9g want %.9g (|diff|=%.3g > %.3g)\n", label,
            got, want, diff, tol);
    ++failures;
  }
}

static void check_true(const char *label, bool cond) {
  if (!cond) {
    fprintf(stderr, "FAIL %s\n", label);
    ++failures;
  }
}

/* Reference scalar from Toeplitz-LS chi2per.c::nll_cybenko_scalar_f. */
static float ref_nll_cybenko_scalar_f(float kappa, const cybenko_model_t *m,
                                      int d) {
  double s = cybenko_kappa_to_S((double)kappa, m);
  if (s <= 0.0 || !double_is_finite_bits(s))
    return 0.0f;

  float k_eff_f = (float)m->k_eff;
  float th_eff_f = (float)m->theta_eff;
  float lg_k_f = (float)m->lgamma_k;
  float s_f = (float)s;

  if (d <= 2) {
    float z = s_f / th_eff_f;
    int order = 2 * d;
    float delta = k_eff_f - (float)order;
    float c = 1.0f;
    float poly = 1.0f;
    float d_poly = 0.0f;
    for (int k = 1; k < order; ++k) {
      c *= (z / (float)k);
      poly += c;
      if (fabsf(delta) > 1.0e-5f) {
        float H_k = 0.0f;
        for (int i = 1; i <= k; ++i)
          H_k += 1.0f / (float)i;
        d_poly += c * (-H_k);
      }
    }
    float nll = z - logf(poly);
    if (fabsf(delta) > 1.0e-5f)
      nll -= delta * (logf(z) + d_poly / poly);
    return (nll > 0.0f && float_is_finite_bits(nll)) ? nll : 0.0f;
  }

  float a = k_eff_f;
  float z = s_f / th_eff_f;
  if (z < a + 3.0f) {
    double c = 1.0;
    double sum = 1.0;
    for (int k = 1; k < 60; ++k) {
      c *= (double)z / ((double)a + (double)(k - 1));
      sum += c;
      if (c < 1.0e-7 * sum)
        break;
    }
    double lg_inc =
        (double)a * log((double)z) - (double)z - log((double)a) + log(sum);
    double nll = (double)lg_k_f - lg_inc;
    return (nll > 0.0 && double_is_finite_bits(nll)) ? (float)nll : 0.0f;
  }

  float inv_z = 1.0f / z;
  float t1 = (a - 1.0f) * inv_z;
  float t2 = t1 * ((a - 2.0f) * inv_z);
  float t3 = t2 * ((a - 3.0f) * inv_z);
  float t4 = t3 * ((a - 4.0f) * inv_z);
  float poly = 1.0f + t1 + t2 + t3 + t4;
  float nll = z - (a - 1.0f) * logf(z) + lg_k_f - logf(poly);
  return (nll > 0.0f && float_is_finite_bits(nll)) ? nll : 0.0f;
}

static void test_cybenko_vector_matches_scalar(void) {
  /* Only d <= 2 is compared against the reference scalar here: the
   * reference scalar's d >= 3 lower-tail branch computes -ln P(a,z),
   * whereas the production vector path computes the periodogram NLL
   * -ln Q(a,z).  The d <= 2 branch is the one changed by this fix. */
  const int degrees[] = {1, 2};
  const int n_effs[] = {100, 1000, 10000};
  const float kappas[] = {1.01f, 1.1f, 2.0f, 10.0f, 100.0f, 1000.0f};

  for (size_t di = 0; di < sizeof(degrees) / sizeof(degrees[0]); ++di) {
    int d = degrees[di];
    for (size_t ni = 0; ni < sizeof(n_effs) / sizeof(n_effs[0]); ++ni) {
      cybenko_model_t m = cybenko_model_init(d, n_effs[ni]);
      for (size_t ki = 0; ki < sizeof(kappas) / sizeof(kappas[0]); ++ki) {
        VEC kv;
        for (int lane = 0; lane < VEC_LEN; ++lane)
          kv.values[lane] = kappas[ki];
        VEC got = nll_cybenko_vec(kv, &m, d);
        float want = ref_nll_cybenko_scalar_f(kappas[ki], &m, d);
        double tol = 2.0e-3 * (1.0 + fabs((double)want));
        char label[128];
        snprintf(label, sizeof(label), "nll_cybenko d=%d n_eff=%d kappa=%.4g",
                 d, n_effs[ni], (double)kappas[ki]);
        check_close(label, (double)got.values[0], (double)want, tol);
      }
    }
  }
}

/* d >= 3 uses the hybrid lower/asymptotic engine; just verify it stays finite
 * and non-negative for representative inputs. */
static void test_cybenko_vector_finite(void) {
  const int degrees[] = {3, 8};
  const int n_effs[] = {100, 1000, 10000};
  const float kappas[] = {1.01f, 1.1f, 2.0f, 10.0f, 100.0f, 1000.0f};
  for (size_t di = 0; di < sizeof(degrees) / sizeof(degrees[0]); ++di) {
    int d = degrees[di];
    for (size_t ni = 0; ni < sizeof(n_effs) / sizeof(n_effs[0]); ++ni) {
      cybenko_model_t m = cybenko_model_init(d, n_effs[ni]);
      for (size_t ki = 0; ki < sizeof(kappas) / sizeof(kappas[0]); ++ki) {
        VEC kv;
        for (int lane = 0; lane < VEC_LEN; ++lane)
          kv.values[lane] = kappas[ki];
        VEC got = nll_cybenko_vec(kv, &m, d);
        for (int lane = 0; lane < VEC_LEN; ++lane) {
          float v = got.values[lane];
          if (!(v >= 0.0f) || !float_is_finite_bits(v)) {
            fprintf(stderr,
                    "FAIL nll_cybenko finite d=%d n_eff=%d kappa=%.4g: %.9g\n",
                    d, n_effs[ni], (double)kappas[ki], (double)v);
            ++failures;
          }
        }
      }
    }
  }
}

/* Build the d=1 GLS sums from a synthetic light curve, then check both the
 * power and the Cybenko condition against the Levinson n=3 path. */
static void test_d1_condition_matches_levinson(void) {
  enum { N = 500 };
  double t[N];
  float y[N], dy[N];
  unsigned seed = 12345u;
  for (int i = 0; i < N; ++i) {
    seed = seed * 1664525u + 1013904223u;
    double u = (double)(seed >> 8) / (double)(1u << 24);
    t[i] = (double)i + 0.1 * u;
    seed = seed * 1664525u + 1013904223u;
    u = (double)(seed >> 8) / (double)(1u << 24);
    y[i] = (float)(sin(2.0 * M_PI * 0.013 * t[i]) + 0.2 * (u - 0.5));
    dy[i] = 0.1f;
  }

  double ws = 0.0, wy = 0.0;
  for (int i = 0; i < N; ++i) {
    double w = 1.0 / ((double)dy[i] * (double)dy[i]);
    ws += w;
    wy += w * (double)y[i];
  }
  double ymean = wy / ws;
  double S = 0.0, C = 0.0, S2 = 0.0, C2 = 0.0, Sh = 0.0, Ch = 0.0,
         chi2_ref = 0.0;
  const double f = 0.013;
  for (int i = 0; i < N; ++i) {
    double w = 1.0 / ((double)dy[i] * (double)dy[i]);
    double yc = (double)y[i] - ymean;
    double ph = 2.0 * M_PI * f * t[i];
    double c1 = cos(ph), s1 = sin(ph);
    double c2 = cos(2.0 * ph), s2 = sin(2.0 * ph);
    C += w * c1;
    S += w * s1;
    C2 += w * c2;
    S2 += w * s2;
    Ch += w * yc * c1;
    Sh += w * yc * s1;
    chi2_ref += w * yc * yc;
  }

  float Sw[3] = {0.0f, (float)S, (float)S2};
  float Cw[3] = {(float)ws, (float)C, (float)C2};
  float Syw[3] = {0.0f, (float)Sh, 0.0f};
  float Cyw[3] = {0.0f, (float)Ch, 0.0f};
  /* yws is sum w*(y-ymean), i.e. zero by construction. */

  float power = 0.0f, cond = 0.0f;
  aov_gls_impl(Sw, Cw, Syw, Cyw, 1, (float)ws, 0.0f, (float)chi2_ref, &power,
               &cond);

  float Rr[3] = {Cw[0], Cw[1], Cw[2]};
  float Ri[3] = {-Sw[0], -Sw[1], -Sw[2]};
  float Yr[3] = {0}, Yi[3] = {0}, Xr[3] = {0}, Xi[3] = {0};
  for (int k = 0; k < 3; ++k) {
    int h = k - 1; /* degree = 1 */
    if (h == 0) {
      Yr[k] = Cyw[0];
      Yi[k] = 0.0f;
    } else if (h > 0) {
      Yr[k] = Cyw[h];
      Yi[k] = -Syw[h];
    } else {
      Yr[k] = Cyw[-h];
      Yi[k] = Syw[-h];
    }
  }
  float Ar[3] = {0}, Ai[3] = {0}, Apr[3] = {0}, Api[3] = {0};
  float lev_cond =
      aov_solve_levinson_scalar(3, Rr, Ri, Yr, Yi, Xr, Xi, Ar, Ai, Apr, Api);
  float lev_power = Xr[0] / (float)chi2_ref;

  check_close("d1 power vs Levinson", (double)power, (double)lev_power, 2.0e-5);
  check_close("d1 condition vs Levinson", (double)cond, (double)lev_cond,
              2.0e-3 * (1.0 + fabs((double)lev_cond)));
}

static void complex_solve(int n, double complex A[9][9], double complex b[9],
                          double complex x[9]) {
  for (int col = 0; col < n; ++col) {
    int pivot = col;
    double best = cabs(A[col][col]);
    for (int row = col + 1; row < n; ++row) {
      double v = cabs(A[row][col]);
      if (v > best) {
        best = v;
        pivot = row;
      }
    }
    if (pivot != col) {
      for (int k = 0; k < n; ++k) {
        double complex tmp = A[col][k];
        A[col][k] = A[pivot][k];
        A[pivot][k] = tmp;
      }
      double complex tb = b[col];
      b[col] = b[pivot];
      b[pivot] = tb;
    }
    double complex d = A[col][col];
    for (int row = col + 1; row < n; ++row) {
      double complex factor = A[row][col] / d;
      for (int k = col; k < n; ++k)
        A[row][k] -= factor * A[col][k];
      b[row] -= factor * b[col];
    }
  }
  for (int row = n - 1; row >= 0; --row) {
    double complex sum = b[row];
    for (int k = row + 1; k < n; ++k)
      sum -= A[row][k] * x[k];
    x[row] = sum / A[row][row];
  }
}

static void test_levinson_quadratic_form(void) {
  for (int d = 1; d <= 4; ++d) {
    int n = 2 * d + 1;
    double complex R[9];
    double mag = 0.0;
    R[0] = 10.0;
    unsigned seed = 9876u + (unsigned)d;
    for (int k = 1; k < n; ++k) {
      seed = seed * 1664525u + 1013904223u;
      double re = ((double)(seed >> 8) / (double)(1u << 24) - 0.5) * 0.2;
      seed = seed * 1664525u + 1013904223u;
      double im = ((double)(seed >> 8) / (double)(1u << 24) - 0.5) * 0.2;
      R[k] = re + im * I;
      mag += cabs(R[k]);
    }
    R[0] = mag + 1.0;

    double complex T[9][9];
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        int lag = i - j;
        T[i][j] = lag >= 0 ? R[lag] : conj(R[-lag]);
      }
    }

    double complex Y[9], X[9];
    for (int i = 0; i < n; ++i) {
      seed = seed * 1664525u + 1013904223u;
      double re = ((double)(seed >> 8) / (double)(1u << 24) - 0.5);
      seed = seed * 1664525u + 1013904223u;
      double im = ((double)(seed >> 8) / (double)(1u << 24) - 0.5);
      Y[i] = re + im * I;
    }
    double complex Acopy[9][9], bcopy[9];
    memcpy(Acopy, T, sizeof(Acopy));
    memcpy(bcopy, Y, sizeof(bcopy));
    complex_solve(n, Acopy, bcopy, X);
    double complex q = 0.0;
    for (int i = 0; i < n; ++i)
      q += conj(Y[i]) * X[i];
    double want = creal(q);

    float Rr[9], Ri[9], Yr[9], Yi[9], Xr[9], Xi[9], Ar[9], Ai[9], Apr[9],
        Api[9];
    for (int i = 0; i < n; ++i) {
      Rr[i] = (float)creal(R[i]);
      Ri[i] = (float)cimag(R[i]);
      Yr[i] = (float)creal(Y[i]);
      Yi[i] = (float)cimag(Y[i]);
    }
    aov_solve_levinson_scalar((size_t)n, Rr, Ri, Yr, Yi, Xr, Xi, Ar, Ai, Apr,
                              Api);
    char label[64];
    snprintf(label, sizeof(label), "Levinson quadratic form d=%d", d);
    check_close(label, (double)Xr[0], want, 2.0e-4 * (1.0 + fabs(want)));
  }
}

static void test_correct_power(void) {
  const float ks[] = {0.01f, 0.1f, 0.5f, 1.0f, 2.0f, 5.0f, 10.0f, 50.0f};
  const float ninvs[] = {1.0e-2f, 1.0e-3f, 1.0e-4f, 1.0e-5f};
  for (size_t ki = 0; ki < sizeof(ks) / sizeof(ks[0]); ++ki) {
    for (size_t ni = 0; ni < sizeof(ninvs) / sizeof(ninvs[0]); ++ni) {
      VEC kv;
      for (int lane = 0; lane < VEC_LEN; ++lane)
        kv.values[lane] = ks[ki];
      VEC got = correctPower(kv, ninvs[ni]);
      float want = correct_power_scalar(ks[ki], ninvs[ni]);
      double tol = 2.0e-3 * (1.0 + fabs((double)want));
      char label[96];
      snprintf(label, sizeof(label), "correctPower K=%.4g nInv=%.4g",
               (double)ks[ki], (double)ninvs[ni]);
      check_close(label, (double)got.values[0], (double)want, tol);
    }
  }

  /* Deliberately invalid bracket: fallback must preserve K. */
  VEC bad;
  for (int lane = 0; lane < VEC_LEN; ++lane)
    bad.values[lane] = 5.0f;
  VEC out = correctPower(bad, 1.0f);
  for (int lane = 0; lane < VEC_LEN; ++lane)
    check_close("correctPower fallback", (double)out.values[lane], 5.0, 1.0e-5);
}

static void test_condition_limit(void) {
  check_true("condition 1e4 not singular", !aov_condition_is_singular(1.0e4f));
  check_true("condition 1e5 singular", aov_condition_is_singular(1.0e5f));

  cybenko_model_t m = cybenko_model_init(1, 1000);
  double at_limit = cybenko_kappa_to_S((double)AOV_CONDITION_LIMIT, &m);
  double above_limit = cybenko_kappa_to_S(1.0e6, &m);
  check_close("kappa clamp uses AOV_CONDITION_LIMIT", above_limit, at_limit,
              1.0e-12);
}

int main(void) {
  test_cybenko_vector_matches_scalar();
  test_cybenko_vector_finite();
  test_d1_condition_matches_levinson();
  test_levinson_quadratic_form();
  test_correct_power();
  test_condition_limit();

  if (failures != 0) {
    fprintf(stderr, "core_fixes_test: %d failure(s)\n", failures);
    return 1;
  }
  printf("core_fixes_test: PASS\n");
  return 0;
}
