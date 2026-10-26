/*
 * Bit-reversal COBRA/COBRAVO and fused-Radix-2 FFT regression/scaling test.
 *
 * Compile once with -DF_BITREV_COBRAVO=0 and once with =1.  Both binaries must
 * print identical permutation and NuFFT hashes; the COBRAVO build must not be
 * slower than COBRA.  The fused Radix-2 pair is checked against the legacy
 * two-pass kernel and against the selected COBRA/COBRAVO permutation.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef NUFFT1_SCALING_BOOTSTRAP
#define NUFFT1_SCALING_BOOTSTRAP 1
#endif
#include "../src/nufft/nufft1.c"

static double now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static float *alloc_float(size_t n) {
  size_t bytes = ((n * sizeof(float) + 63U) / 64U) * 64U;
  return (float *)aligned_alloc(64, bytes);
}

static int check_order(uint32_t order) {
  size_t n = (size_t)1U << order;
  float *ref_r = alloc_float(n), *ref_i = alloc_float(n);
  float *got_r = alloc_float(n), *got_i = alloc_float(n);
  float *a_r = alloc_float(2U * (size_t)BLOCK_ELEMS);
  float *a_i = alloc_float(2U * (size_t)BLOCK_ELEMS);
  float *b_r = alloc_float((size_t)BLOCK_ELEMS);
  float *b_i = alloc_float((size_t)BLOCK_ELEMS);
  if (!ref_r || !ref_i || !got_r || !got_i || !a_r || !a_i || !b_r || !b_i) {
    fprintf(stderr, "alloc failed at N=%zu\n", n);
    return 1;
  }
  for (size_t k = 0; k < n; ++k) {
    ref_r[k] = got_r[k] = (float)k + 0.25f;
    ref_i[k] = got_i[k] = -(float)k - 0.5f;
  }
  table_shuffle(ref_r, ref_i, order);

  cobra_shuffle(got_r, got_i, order, a_r, a_i);
  int cobra_ok = memcmp(ref_r, got_r, n * sizeof(float)) == 0 &&
                 memcmp(ref_i, got_i, n * sizeof(float)) == 0;

  for (size_t k = 0; k < n; ++k) {
    got_r[k] = (float)k + 0.25f;
    got_i[k] = -(float)k - 0.5f;
  }
  cobravo_init();
  cobravo_shuffle(got_r, got_i, order, a_r, a_i, b_r, b_i);
  int cobravo_ok = memcmp(ref_r, got_r, n * sizeof(float)) == 0 &&
                   memcmp(ref_i, got_i, n * sizeof(float)) == 0;

  for (size_t k = 0; k < n; ++k) {
    got_r[k] = (float)k + 0.25f;
    got_i[k] = -(float)k - 0.5f;
  }
  bit_reverse_permutation(got_r, got_i, (uint32_t)n, a_r, a_i);
  int selected_ok = memcmp(ref_r, got_r, n * sizeof(float)) == 0 &&
                    memcmp(ref_i, got_i, n * sizeof(float)) == 0;

  printf("CHECK order=%u N=%zu cobra=%s cobravo=%s selected=%s\n", order, n,
         cobra_ok ? "OK" : "FAIL", cobravo_ok ? "OK" : "FAIL",
         selected_ok ? "OK" : "FAIL");
  free(ref_r);
  free(ref_i);
  free(got_r);
  free(got_i);
  free(a_r);
  free(a_i);
  free(b_r);
  free(b_i);
  return (cobra_ok && cobravo_ok && selected_ok) ? 0 : 1;
}

static double perm_time(int use_cobravo, uint32_t order, int reps) {
  size_t n = (size_t)1U << order;
  float *r = alloc_float(n), *i = alloc_float(n);
  float *a_r = alloc_float(2U * (size_t)BLOCK_ELEMS);
  float *a_i = alloc_float(2U * (size_t)BLOCK_ELEMS);
  float *b_r = alloc_float((size_t)BLOCK_ELEMS);
  float *b_i = alloc_float((size_t)BLOCK_ELEMS);
  if (!r || !i || !a_r || !a_i || !b_r || !b_i)
    return NAN;
  for (size_t k = 0; k < n; ++k) {
    r[k] = (float)k;
    i[k] = -(float)k;
  }
  cobravo_init();
  for (int w = 0; w < 3; ++w) {
    if (use_cobravo)
      cobravo_shuffle(r, i, order, a_r, a_i, b_r, b_i);
    else
      cobra_shuffle(r, i, order, a_r, a_i);
  }
  double best = 1e30;
  for (int rep = 0; rep < reps; ++rep) {
    double t0 = now_seconds();
    if (use_cobravo)
      cobravo_shuffle(r, i, order, a_r, a_i, b_r, b_i);
    else
      cobra_shuffle(r, i, order, a_r, a_i);
    double dt = now_seconds() - t0;
    if (dt < best)
      best = dt;
  }
  free(r);
  free(i);
  free(a_r);
  free(a_i);
  free(b_r);
  free(b_i);
  return best;
}

static uint64_t hash_floats(const float *r, const float *i, size_t n) {
  uint64_t h = 1469598103934665603ULL;
  for (size_t k = 0; k < n; ++k) {
    uint32_t u;
    memcpy(&u, &r[k], sizeof(u));
    h = (h ^ u) * 1099511628211ULL;
    memcpy(&u, &i[k], sizeof(u));
    h = (h ^ u) * 1099511628211ULL;
  }
  return h;
}

static int pipeline_case(int N, int M) {
  nufft1_external_sizes sizes =
      nufft1_external_sizes_for_plan(N, NUFFT1_PSWF43);
  float *tw_r = alloc_float(sizes.twiddle_len),
        *tw_i = alloc_float(sizes.twiddle_len);
  float *fft_r = alloc_float(sizes.fft_len),
        *fft_i = alloc_float(sizes.fft_len);
  float *cob_r = alloc_float(sizes.cobra_len),
        *cob_i = alloc_float(sizes.cobra_len);
  nufft1_plan *plan =
      nufft1_initialize_shared(M, N, 1, NUFFT1_PSWF43, tw_r, tw_i);
  nufft1_workspace *ws =
      nufft1_create_workspace(plan, fft_r, fft_i, cob_r, cob_i);
  if (!plan || !ws) {
    nufft1_free_workspace(ws);
    nufft1_free_plan(plan);
    free(tw_r);
    free(tw_i);
    free(fft_r);
    free(fft_i);
    free(cob_r);
    free(cob_i);
    return 1;
  }

  double *x = (double *)malloc((size_t)M * sizeof(double));
  float *yr = (float *)malloc((size_t)M * sizeof(float));
  float *yi = (float *)malloc((size_t)M * sizeof(float));
  srand(1);
  for (int k = 0; k < M; ++k) {
    x[k] = (double)rand() / (double)RAND_MAX;
    yr[k] = (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
    yi[k] = (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
  }
  nufft1_precompute(ws, x, M, 1.0);
  int nout = nufft1_output_len(plan);
  float *out_r = (float *)calloc((size_t)nout, sizeof(float));
  float *out_i = (float *)calloc((size_t)nout, sizeof(float));
  nufft1_execute(ws, yr, yi, out_r, out_i, 1);
  double best = 1e30;
  for (int rep = 0; rep < 10; ++rep) {
    double t0 = now_seconds();
    nufft1_execute(ws, yr, yi, out_r, out_i, 1);
    double dt = now_seconds() - t0;
    if (dt < best)
      best = dt;
  }
  printf("HASH N=%d Nfft=%d M=%d hash=%016llx exec_us=%.2f\n", N,
         nufft1_plan_len(plan) * 2, M,
         (unsigned long long)hash_floats(out_r, out_i, (size_t)nout),
         best * 1e6);
  free(x);
  free(yr);
  free(yi);
  free(out_r);
  free(out_i);
  nufft1_free_workspace(ws);
  nufft1_free_plan(plan);
  free(tw_r);
  free(tw_i);
  free(fft_r);
  free(fft_i);
  free(cob_r);
  free(cob_i);
  return 0;
}

static int check_fused_fft(uint32_t order) {
  size_t n = (size_t)1U << order;
  float *ref_r = alloc_float(n), *ref_i = alloc_float(n);
  float *got_r = alloc_float(n), *got_i = alloc_float(n);
  float *tw_old_r = alloc_float(n), *tw_old_i = alloc_float(n);
  float *tw_fused_r = alloc_float(n), *tw_fused_i = alloc_float(n);
  if (!ref_r || !ref_i || !got_r || !got_i || !tw_old_r || !tw_old_i ||
      !tw_fused_r || !tw_fused_i) {
    fprintf(stderr, "alloc failed at N=%zu\n", n);
    return 1;
  }
  srand(900U + order);
  for (size_t k = 0; k < n; ++k) {
    ref_r[k] = got_r[k] =
        (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
    ref_i[k] = got_i[k] =
        (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
    tw_old_r[k] = tw_old_i[k] = tw_fused_r[k] = tw_fused_i[k] = 0.0f;
  }
  generate_fft_buffer((uint32_t)n, tw_old_r, tw_old_i);
  generate_fft_buffer_fused_r2((uint32_t)n, tw_fused_r, tw_fused_i);
  sande_tukey_in_place(ref_r, ref_i, tw_old_r, tw_old_i, (uint32_t)n);
  sande_tukey_in_place_fused_r2(got_r, got_i, tw_fused_r, tw_fused_i,
                                (uint32_t)n);
  table_shuffle(ref_r, ref_i, order);
  table_shuffle(got_r, got_i, order);

  double max_abs = 0.0, max_rel = 0.0;
  for (size_t k = 0; k < n; ++k) {
    double dr = fabs((double)got_r[k] - (double)ref_r[k]);
    double di = fabs((double)got_i[k] - (double)ref_i[k]);
    double mag = sqrt((double)ref_r[k] * (double)ref_r[k] +
                      (double)ref_i[k] * (double)ref_i[k]);
    if (dr > max_abs)
      max_abs = dr;
    if (di > max_abs)
      max_abs = di;
    if (mag > 1e-9) {
      double rel = sqrt(dr * dr + di * di) / mag;
      if (rel > max_rel)
        max_rel = rel;
    }
  }
  int ok = (max_abs < 5e-3 && max_rel < 5e-4);
  printf("FUSED order=%u N=%zu max_rel=%.3e %s\n", order, n, max_rel,
         ok ? "OK" : "FAIL");
  free(ref_r);
  free(ref_i);
  free(got_r);
  free(got_i);
  free(tw_old_r);
  free(tw_old_i);
  free(tw_fused_r);
  free(tw_fused_i);
  return ok ? 0 : 1;
}

static int check_fused_perm(uint32_t order) {
  size_t n = (size_t)1U << order;
  float *ref_r = alloc_float(n), *ref_i = alloc_float(n);
  float *got_r = alloc_float(n), *got_i = alloc_float(n);
  float *tw_old_r = alloc_float(n), *tw_old_i = alloc_float(n);
  float *tw_fused_r = alloc_float(n), *tw_fused_i = alloc_float(n);
  float *buf_r = alloc_float(2U * (size_t)BLOCK_ELEMS);
  float *buf_i = alloc_float(2U * (size_t)BLOCK_ELEMS);
  if (!ref_r || !ref_i || !got_r || !got_i || !tw_old_r || !tw_old_i ||
      !tw_fused_r || !tw_fused_i || !buf_r || !buf_i) {
    fprintf(stderr, "alloc failed at N=%zu\n", n);
    return 1;
  }
  srand(1800U + order);
  for (size_t k = 0; k < n; ++k) {
    ref_r[k] = got_r[k] =
        (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
    ref_i[k] = got_i[k] =
        (float)(2.0 * (double)rand() / (double)RAND_MAX - 1.0);
  }
  generate_fft_buffer((uint32_t)n, tw_old_r, tw_old_i);
  generate_fft_buffer_fused_r2((uint32_t)n, tw_fused_r, tw_fused_i);
  sande_tukey_in_place(ref_r, ref_i, tw_old_r, tw_old_i, (uint32_t)n);
  sande_tukey_in_place_fused_r2(got_r, got_i, tw_fused_r, tw_fused_i,
                                (uint32_t)n);
  table_shuffle(ref_r, ref_i, order);
  bit_reverse_permutation(got_r, got_i, (uint32_t)n, buf_r, buf_i);

  double max_abs = 0.0, max_rel = 0.0;
  for (size_t k = 0; k < n; ++k) {
    double dr = fabs((double)got_r[k] - (double)ref_r[k]);
    double di = fabs((double)got_i[k] - (double)ref_i[k]);
    double mag = sqrt((double)ref_r[k] * (double)ref_r[k] +
                      (double)ref_i[k] * (double)ref_i[k]);
    if (dr > max_abs)
      max_abs = dr;
    if (di > max_abs)
      max_abs = di;
    if (mag > 1e-9) {
      double rel = sqrt(dr * dr + di * di) / mag;
      if (rel > max_rel)
        max_rel = rel;
    }
  }
  int ok = (max_abs < 5e-3 && max_rel < 5e-4);
  printf("FUSEDPERM order=%u N=%zu backend=%s max_rel=%.3e %s\n", order, n,
         F_BITREV_COBRAVO ? "COBRAVO" : "COBRA", max_rel, ok ? "OK" : "FAIL");
  free(ref_r);
  free(ref_i);
  free(got_r);
  free(got_i);
  free(tw_old_r);
  free(tw_old_i);
  free(tw_fused_r);
  free(tw_fused_i);
  free(buf_r);
  free(buf_i);
  return ok ? 0 : 1;
}

int main(void) {
  printf("F_BITREV_COBRAVO=%d gather_lanes=%d\n", F_BITREV_COBRAVO,
         COBRAVO_GATHER_LANES);
  int failed = 0;
  for (uint32_t order = 13; order <= 18; ++order)
    failed |= check_order(order);
  for (uint32_t order = 13; order <= 18; ++order)
    failed |= check_fused_fft(order);
  for (uint32_t order = 13; order <= 18; ++order)
    failed |= check_fused_perm(order);

  double t13_c = perm_time(0, 13, 40), t13_v = perm_time(1, 13, 40);
  double t16_c = perm_time(0, 16, 25), t16_v = perm_time(1, 16, 25);
  printf("PERM N=8192 cobra_us=%.3f cobravo_us=%.3f speedup=%.2f\n",
         t13_c * 1e6, t13_v * 1e6, t13_c / t13_v);
  printf("PERM N=65536 cobra_us=%.3f cobravo_us=%.3f speedup=%.2f\n",
         t16_c * 1e6, t16_v * 1e6, t16_c / t16_v);
  if (!(t16_v < t16_c * 1.05)) {
    fprintf(stderr, "COBRAVO is not faster at N=65536\n");
    failed = 1;
  }

  failed |= pipeline_case(4096, 128);
  failed |= pipeline_case(8192, 128);
  failed |= pipeline_case(16384, 128);
  failed |= pipeline_case(32768, 128);
  failed |= pipeline_case(65536, 128);
  return failed ? 1 : 0;
}
