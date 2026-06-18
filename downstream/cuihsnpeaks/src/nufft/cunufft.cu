#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cunufft.h"

/* -------------------------------------------------------------------------
 * Trigonometric Remez Approximations (matching trig.h)
 * ------------------------------------------------------------------------- */

__host__ __device__ static inline float sin2pif_scalar(float x) {
    float f = x - (float)((int)x);
    if (f < 0.0f) f += 1.0f;
    float sign = 1.0f;
    if (f >= 0.5f) {
        sign = -1.0f;
        f -= 0.5f;
    }
    if (f > 0.25f) f = 0.5f - f;
    float f2 = f * f;
    float p = f2 * 39.536706065730207835108712734262f - 76.549782293595742666226937116116f;
    p = p * f2 + 81.601004073261773523492199897936f;
    p = p * f2 - 41.341655031416278077153126232486f;
    p = p * f2 + 6.2831851600894774430188071795666f;
    p *= f;
    return p * sign;
}

__host__ __device__ static inline float cos2pif_scalar(float x) {
    float f = x - (float)((int)x);
    if (f < 0.0f) f += 1.0f;
    if (f > 0.5f) f = 1.0f - f;
    float sign = 1.0f;
    if (f > 0.25f) {
        sign = -1.0f;
        f = 0.5f - f;
    }
    float f2 = f * f;
    float p = f2 * 56.242380464873243259663276802701f - 85.240330322699427859509454517828f;
    p = p * f2 + 64.934590626780991246193352727536f;
    p = p * f2 - 19.739171434702393618770795066531f;
    p = p * f2 + 0.99999995346667013630639784578184f;
    return p * sign;
}

/* -------------------------------------------------------------------------
 * Remez Exponential Approximation (matching nufft1.c)
 * ------------------------------------------------------------------------- */

__host__ __device__ static inline float expf_remez_scalar(float x) {
    float y = x * 1.4426950408889634f;
    int i = (int)y;
    float f = y - (float)i;
    if (f < 0.0f) { i -= 1; f += 1.0f; }
    float p = 0.0018964611454333f;
    p = p * f + 0.0089428289841091f;
    p = p * f + 0.0558662463045207f;
    p = p * f + 0.2401397110907695f;
    p = p * f + 0.6931547524751674f;
    p = p * f + 0.9999998931108267f;
    union { int i; float f; } bc;
    bc.i = (i + 127) << 23;
    return p * bc.f;
}

/* -------------------------------------------------------------------------
 * PSWF Evaluation (matching nufft1.c)
 * ------------------------------------------------------------------------- */

__host__ __device__ static inline float pswf0_scalar(float z, cunufft_mode mode) {
    float c = (mode == CUNUFFT_PSWF43) ? 15.658f
              : (float)(3.14159265358979323846 * 8.0 * 0.75 - 0.05);
    float z2 = z * z;
    float poly;
    if (mode == CUNUFFT_PSWF43) {
        poly = 6.059434412e-01f;
        poly = fmaf(poly, z2, -3.077828580e+00f);
        poly = fmaf(poly, z2,  3.170650492e+00f);
        poly = fmaf(poly, z2,  3.081707523e-01f);
        poly = fmaf(poly, z2,  8.972535502e-01f);
        poly = fmaf(poly, z2, -1.574809829e+00f);
        poly = fmaf(poly, z2, -1.705772933e+00f);
        poly = fmaf(poly, z2,  3.814843807e-01f);
        poly = fmaf(poly, z2,  1.000000159e+00f);
    } else {
        poly = 2.1704596899999999f;
        poly = fmaf(poly, z2, -7.1979701699999996f);
        poly = fmaf(poly, z2,  5.6755215970000004f);
        poly = fmaf(poly, z2,  0.41027796129999999f);
        poly = fmaf(poly, z2,  1.5911306080000001f);
        poly = fmaf(poly, z2, -1.930650019f);
        poly = fmaf(poly, z2, -2.099454572f);
        poly = fmaf(poly, z2,  0.38030810120000003f);
        poly = fmaf(poly, z2,  1.0000001350000001f);
    }
    return expf_remez_scalar(-0.5f * c * z2) * poly;
}

/* -------------------------------------------------------------------------
 * Structures
 * ------------------------------------------------------------------------- */

struct cunufft_plan {
    int Mpoints;
    int N;
    int Nout;
    int Nfft;
    int spread_stride;
    int fft_alloc_len;
    int output_shift;
    int num_factors;
    cunufft_mode mode;
    float alpha;
    float *d_deconv;
    float *d_twiddle_real;
    float *d_twiddle_imag;
};

struct cunufft_workspace {
    const cunufft_plan *plan;
    int active_mpoints;
    /* Precomputed data for ALL factors */
    int   *d_spread_base_idx;  /* [num_factors * Mpoints]                    */
    float *d_spread_weight;    /* [num_factors * Mpoints * spread_stride]    */
    float *d_shift_r;          /* [num_factors * Mpoints]                    */
    float *d_shift_i;          /* [num_factors * Mpoints]                    */
    /* Single FFT working buffer (reused per factor) */
    float *d_fft_real;         /* [fft_alloc_len]                            */
    float *d_fft_imag;         /* [fft_alloc_len]                            */
    /* Shared input */
    float *d_y_real;           /* [Mpoints]                                  */
    float *d_y_imag;           /* [Mpoints]                                  */
    /* Batched output */
    float *d_out_real;         /* [num_factors * Nout]                       */
    float *d_out_imag;         /* [num_factors * Nout]                       */
};

/* -------------------------------------------------------------------------
 * Host helpers
 * ------------------------------------------------------------------------- */

static inline int bitceil_int(unsigned int n) {
    if (n <= 1) return 1;
    unsigned int v = n - 1U;
    int shift = 0;
    while (v) { v >>= 1; ++shift; }
    return 1 << shift;
}

static inline uint32_t intlog2(uint32_t input) {
    uint32_t output;
    frexp(input >> 1, (int *)&output);
    return output;
}

static void get_lege_roots(int n, float *roots, float *weights) {
    int half_n = (n + 1) / 2;
    for (int i = 1; i <= half_n; ++i) {
        float x = cos2pif_scalar(0.5f * (float)(i - 0.25) / (float)(n + 0.5));
        float p1 = 0.0f, p2 = 0.0f, dp = 0.0f, step;
        int iter = 0;
        do {
            p1 = 1.0f; p2 = 0.0f;
            for (int j = 1; j <= n; ++j) {
                float p0 = (((float)(2*j-1)*x*p1) - ((float)(j-1)*p2)) / (float)j;
                p2 = p1; p1 = p0;
            }
            dp = (float)n * (x*p1 - p2) / (x*x - 1.0f);
            step = p1 / dp;
            x -= step;
        } while (fabsf(step) > 1.0e-7f && ++iter < 64);
        roots[i-1]   = x;  weights[i-1]   = 2.0f / ((1.0f - x*x) * dp * dp);
        roots[n-i]   = -x; weights[n-i]   = weights[i-1];
    }
}

static int initialize_deconv(cunufft_plan *plan, float *h_deconv) {
    int p = plan->mode == CUNUFFT_PSWF43 ? (4*8+16) : (int)(1.5*8+2);
    int num_gl = 2 * p;
    float *gl_n = (float *)malloc(num_gl * sizeof(float));
    float *gl_w = (float *)malloc(num_gl * sizeof(float));
    float *pv   = (float *)malloc(num_gl * sizeof(float));
    if (!gl_n || !gl_w || !pv) { free(pv); free(gl_n); free(gl_w); return -3; }
    get_lege_roots(num_gl, gl_n, gl_w);
    for (int i = 0; i < p; ++i) pv[i] = pswf0_scalar(gl_n[i], plan->mode);
    for (int j = 0; j < p; ++j) pv[j] *= gl_w[j];
    for (int k = 0; k < plan->Nout; ++k) {
        float xi_k = plan->alpha * (float)(k - plan->output_shift);
        float phi = 0.0f;
        for (int j = 0; j < p; ++j) phi += pv[j] * cos2pif_scalar(xi_k * gl_n[j]);
        h_deconv[k] = 1.0f / (8.0f * phi);
    }
    free(pv); free(gl_n); free(gl_w);
    return 0;
}

static void generate_fft_buffer(uint32_t n, float *rb, float *ib) {
    uint32_t shift = 0;
    for (uint32_t step = n; step > 1U; step >>= 1U) {
        uint32_t hs = step >> 1U;
        for (uint32_t j = 0; j < hs; ++j) {
            float angle = -(float)j / (float)step;
            rb[shift+j] = cos2pif_scalar(angle);
            ib[shift+j] = -sin2pif_scalar(angle);
        }
        shift += hs;
    }
}

/* -------------------------------------------------------------------------
 * Device helpers
 * ------------------------------------------------------------------------- */

__device__ static inline int positive_mod_gpu(int value, int modulus) {
    int result = value % modulus;
    return result < 0 ? result + modulus : result;
}

/* -------------------------------------------------------------------------
 * CUDA Kernels (single-factor, reusable per iteration)
 * ------------------------------------------------------------------------- */

__global__ void cunufft_precompute_kernel(
    const double * __restrict__ d_x, int Mpoints, double df, int num_factors,
    int Nfft, int output_shift, int spread_stride, cunufft_mode mode,
    int *d_spread_base_idx, float *d_spread_weight,
    float *d_shift_r, float *d_shift_i)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= Mpoints * num_factors) return;
    int f_idx = idx / Mpoints;
    int m     = idx % Mpoints;
    int k_factor = f_idx + 1;
    double xm = d_x[m] * df * (double)k_factor;
    double phase = xm * (double)output_shift;
    double phase_frac = phase - (double)((int)phase);
    d_shift_r[idx] = cos2pif_scalar((float)phase_frac);
    d_shift_i[idx] = sin2pif_scalar((float)phase_frac);
    double x_n = xm * (double)Nfft;
    int m_left = (int)ceil(x_n - 4.0);
    d_spread_base_idx[idx] = positive_mod_gpu(m_left, Nfft);
    size_t w_base = (size_t)idx * (size_t)spread_stride;
    for (int l = 0; l < spread_stride; ++l) {
        double dist = (double)(m_left + l) - x_n;
        d_spread_weight[w_base + l] = pswf0_scalar((float)(2.0 * dist / 8.0), mode);
    }
}

/* Spread: one factor at a time, parallelised over (point * spread_stride) */
__global__ void cunufft_spread_kernel(
    const float * __restrict__ d_y_real, const float * __restrict__ d_y_imag,
    const float * __restrict__ d_shift_r, const float * __restrict__ d_shift_i,
    const int   * __restrict__ d_spread_base_idx,
    const float * __restrict__ d_spread_weight,
    int active_mpoints, int spread_stride, int Nfft,
    float *d_fft_real, float *d_fft_imag)
{
    int tid   = blockIdx.x * blockDim.x + threadIdx.x;
    int total = active_mpoints * spread_stride;
    if (tid >= total) return;
    int m = tid >> 3;
    int l = tid & 7;
    float yr = d_y_real[m], yi = d_y_imag[m];
    float sr = d_shift_r[m], si = d_shift_i[m];
    float sy_r = yr * sr - yi * si;
    float sy_i = yr * si + yi * sr;
    int base = d_spread_base_idx[m];
    float kval = d_spread_weight[((size_t)m << 3) + l];
    int idx = base + l;
    if (idx >= Nfft) idx -= Nfft;
    atomicAdd(&d_fft_real[idx], sy_r * kval);
    atomicAdd(&d_fft_imag[idx], sy_i * kval);
}

/* Sande-Tukey butterfly: global memory, bitwise indexing */
__global__ void sande_tukey_stage_kernel(
    float *real, float *imag,
    const float * __restrict__ tw_r, const float * __restrict__ tw_i,
    uint32_t half_step, uint32_t log_half_step, uint32_t shift, uint32_t n)
{
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= (n >> 1)) return;
    uint32_t i = (tid >> log_half_step) << (log_half_step + 1);
    uint32_t j = tid & (half_step - 1);
    uint32_t ie = i + j, io = ie + half_step;
    float re = real[ie], ime = imag[ie], ro = real[io], imo = imag[io];
    real[ie] = re + ro;  imag[ie] = ime + imo;
    float rt = re - ro, it = ime - imo;
    float tr = tw_r[shift+j], ti = tw_i[shift+j];
    real[io] = rt*tr - it*ti;  imag[io] = rt*ti + it*tr;
}

/* Sande-Tukey: shared-memory stages for step <= 1024 */
__global__ void sande_tukey_shared_kernel(
    float *real, float *imag,
    const float * __restrict__ tw_r, const float * __restrict__ tw_i,
    uint32_t start_step, uint32_t log_start_step, uint32_t start_shift, uint32_t n)
{
    extern __shared__ float s[];
    float *sr = s, *si = s + start_step;
    uint32_t bo = blockIdx.x * start_step, tid = threadIdx.x;
    uint32_t hs0 = start_step >> 1;
    sr[tid]       = real[bo + tid];
    si[tid]       = imag[bo + tid];
    sr[tid + hs0] = real[bo + tid + hs0];
    si[tid + hs0] = imag[bo + tid + hs0];
    __syncthreads();
    uint32_t shift = start_shift, lhs = log_start_step - 1;
    for (uint32_t step = start_step; step > 1U; step >>= 1U) {
        uint32_t hs = step >> 1U;
        uint32_t ii = (tid >> lhs) << (lhs + 1);
        uint32_t jj = tid & (hs - 1);
        uint32_t se = ii + jj, so = se + hs;
        float re = sr[se], ie = si[se], ro = sr[so], io = si[so];
        sr[se] = re + ro;  si[se] = ie + io;
        float rt = re - ro, it = ie - io;
        float tr = tw_r[shift+jj], ti = tw_i[shift+jj];
        sr[so] = rt*tr - it*ti;  si[so] = rt*ti + it*tr;
        shift += hs; lhs--;
        __syncthreads();
    }
    real[bo + tid]       = sr[tid];
    imag[bo + tid]       = si[tid];
    real[bo + tid + hs0] = sr[tid + hs0];
    imag[bo + tid + hs0] = si[tid + hs0];
}

/* Bit-reverse using hardware __brev */
__global__ void bit_reverse_kernel(float *real, float *imag, uint32_t n, uint32_t log_n) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    uint32_t j = __brev(tid) >> (32 - log_n);
    if (j > tid) {
        float tr = real[tid], ti = imag[tid];
        real[tid] = real[j];  imag[tid] = imag[j];
        real[j] = tr;         imag[j] = ti;
    }
}

/* Deconvolution scaling: reads single FFT buffer, writes to output slice */
__global__ void cunufft_deconv_kernel(
    const float * __restrict__ fft_r, const float * __restrict__ fft_i,
    const float * __restrict__ deconv, int Nfft, int Nout, int output_shift,
    float *out_r, float *out_i)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= Nout) return;
    int kp = k - output_shift;
    int fi = kp >= 0 ? kp : Nfft + kp;
    float s = deconv[k];
    out_r[k] = fft_r[fi] * s;
    out_i[k] = fft_i[fi] * s;
}

/* -------------------------------------------------------------------------
 * Public API
 * ------------------------------------------------------------------------- */

cunufft_plan *cunufft_initialize(int Mpoints, int N, int freq_factor, cunufft_mode mode) {
    if (Mpoints <= 0 || N <= 0 || freq_factor <= 0) return NULL;
    if (mode != CUNUFFT_PSWF21 && mode != CUNUFFT_PSWF43) return NULL;
    if (mode == CUNUFFT_PSWF43 && (N % 4) != 0) return NULL;
    cunufft_plan *p = (cunufft_plan *)calloc(1, sizeof(*p));
    if (!p) return NULL;
    p->Mpoints      = Mpoints;
    p->N             = N;
    p->Nout          = mode == CUNUFFT_PSWF43 ? N + (N >> 1) : N;
    p->Nfft          = bitceil_int(2 * N);
    p->output_shift  = mode == CUNUFFT_PSWF43 ? 3*N/4 : N/2;
    p->spread_stride = 8;
    p->fft_alloc_len = p->Nfft;
    p->num_factors   = freq_factor;
    p->mode          = mode;
    p->alpha         = 4.0f / (float)p->Nfft;
    float *hd  = (float *)malloc(p->Nout * sizeof(float));
    float *htr = (float *)malloc(p->Nfft * sizeof(float));
    float *hti = (float *)malloc(p->Nfft * sizeof(float));
    if (!hd || !htr || !hti) { free(hd); free(htr); free(hti); free(p); return NULL; }
    initialize_deconv(p, hd);
    generate_fft_buffer((uint32_t)p->Nfft, htr, hti);
    cudaMalloc(&p->d_deconv,       p->Nout * sizeof(float));
    cudaMalloc(&p->d_twiddle_real,  p->Nfft * sizeof(float));
    cudaMalloc(&p->d_twiddle_imag,  p->Nfft * sizeof(float));
    cudaMemcpy(p->d_deconv,       hd,  p->Nout * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(p->d_twiddle_real, htr, p->Nfft * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(p->d_twiddle_imag, hti, p->Nfft * sizeof(float), cudaMemcpyHostToDevice);
    free(hd); free(htr); free(hti);
    return p;
}

cunufft_workspace *cunufft_create_workspace(const cunufft_plan *plan) {
    if (!plan) return NULL;
    cunufft_workspace *ws = (cunufft_workspace *)calloc(1, sizeof(*ws));
    if (!ws) return NULL;
    ws->plan = plan;
    int ne = plan->Mpoints * plan->num_factors;
    cudaMalloc(&ws->d_spread_base_idx, ne * sizeof(int));
    cudaMalloc(&ws->d_spread_weight,   ne * plan->spread_stride * sizeof(float));
    cudaMalloc(&ws->d_shift_r,         ne * sizeof(float));
    cudaMalloc(&ws->d_shift_i,         ne * sizeof(float));
    /* Single FFT buffer — reused each iteration */
    cudaMalloc(&ws->d_fft_real, plan->fft_alloc_len * sizeof(float));
    cudaMalloc(&ws->d_fft_imag, plan->fft_alloc_len * sizeof(float));
    cudaMalloc(&ws->d_y_real, plan->Mpoints * sizeof(float));
    cudaMalloc(&ws->d_y_imag, plan->Mpoints * sizeof(float));
    size_t out_n = (size_t)plan->num_factors * (size_t)plan->Nout;
    cudaMalloc(&ws->d_out_real, out_n * sizeof(float));
    cudaMalloc(&ws->d_out_imag, out_n * sizeof(float));
    return ws;
}

void cunufft_precompute(cunufft_workspace *ws, const double *x, int Mpoints, double df) {
    if (!ws || !ws->plan || !x || Mpoints <= 0) return;
    const cunufft_plan *p = ws->plan;
    ws->active_mpoints = Mpoints;
    double *dx; cudaMalloc(&dx, Mpoints * sizeof(double));
    cudaMemcpy(dx, x, Mpoints * sizeof(double), cudaMemcpyHostToDevice);
    int total = Mpoints * p->num_factors;
    cunufft_precompute_kernel<<<(total+255)/256, 256>>>(
        dx, Mpoints, df, p->num_factors, p->Nfft, p->output_shift,
        p->spread_stride, p->mode,
        ws->d_spread_base_idx, ws->d_spread_weight, ws->d_shift_r, ws->d_shift_i);
    cudaDeviceSynchronize();
    cudaFree(dx);
}

float cunufft_execute_batch(const cunufft_workspace *ws,
                            const float *y_real, const float *y_imag,
                            float *out_real, float *out_imag)
{
    if (!ws || !ws->plan || !y_real || !y_imag) return 0.0f;
    const cunufft_plan *p = ws->plan;
    if (ws->active_mpoints <= 0) return 0.0f;
    int M  = ws->active_mpoints;
    int ss = p->spread_stride;

    /* Upload shared input */
    cudaMemcpy(ws->d_y_real, y_real, M * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ws->d_y_imag, y_imag, M * sizeof(float), cudaMemcpyHostToDevice);

    cudaEvent_t t0, t1;
    cudaEventCreate(&t0); cudaEventCreate(&t1);
    cudaEventRecord(t0);

    uint32_t nfft  = (uint32_t)p->Nfft;
    uint32_t log_n = intlog2(nfft);

    for (int fi = 0; fi < p->num_factors; ++fi) {
        size_t off   = (size_t)fi * (size_t)p->Mpoints;
        size_t off_w = off * (size_t)ss;

        /* 1. Zero FFT grid */
        cudaMemsetAsync(ws->d_fft_real, 0, nfft * sizeof(float));
        cudaMemsetAsync(ws->d_fft_imag, 0, nfft * sizeof(float));

        /* 2. Spread */
        int total_sp = M * ss;
        cunufft_spread_kernel<<<(total_sp+255)/256, 256>>>(
            ws->d_y_real, ws->d_y_imag,
            ws->d_shift_r + off, ws->d_shift_i + off,
            ws->d_spread_base_idx + off, ws->d_spread_weight + off_w,
            M, ss, (int)nfft, ws->d_fft_real, ws->d_fft_imag);

        /* 3. Sande-Tukey FFT */
        uint32_t shift = 0, step = nfft;
        while (step > 1024) {
            uint32_t hs  = step >> 1;
            uint32_t lhs = intlog2(hs);
            sande_tukey_stage_kernel<<<((nfft/2)+255)/256, 256>>>(
                ws->d_fft_real, ws->d_fft_imag,
                p->d_twiddle_real, p->d_twiddle_imag,
                hs, lhs, shift, nfft);
            shift += hs; step = hs;
        }
        if (step > 1) {
            uint32_t nb = nfft / step, tpb = step / 2;
            sande_tukey_shared_kernel<<<nb, tpb, 2*step*sizeof(float)>>>(
                ws->d_fft_real, ws->d_fft_imag,
                p->d_twiddle_real, p->d_twiddle_imag,
                step, intlog2(step), shift, nfft);
        }

        /* 4. Bit-reverse */
        bit_reverse_kernel<<<(nfft+255)/256, 256>>>(
            ws->d_fft_real, ws->d_fft_imag, nfft, log_n);

        /* 5. Deconv -> output slice */
        cunufft_deconv_kernel<<<(p->Nout+255)/256, 256>>>(
            ws->d_fft_real, ws->d_fft_imag, p->d_deconv,
            p->Nfft, p->Nout, p->output_shift,
            ws->d_out_real + (size_t)fi * (size_t)p->Nout,
            ws->d_out_imag + (size_t)fi * (size_t)p->Nout);
    }

    cudaEventRecord(t1);
    cudaEventSynchronize(t1);
    float ms = 0.0f;
    cudaEventElapsedTime(&ms, t0, t1);

    size_t out_n = (size_t)p->num_factors * (size_t)p->Nout;
    cudaMemcpy(out_real, ws->d_out_real, out_n * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(out_imag, ws->d_out_imag, out_n * sizeof(float), cudaMemcpyDeviceToHost);
    cudaEventDestroy(t0); cudaEventDestroy(t1);
    return ms;
}

void cunufft_free_workspace(cunufft_workspace *ws) {
    if (!ws) return;
    cudaFree(ws->d_spread_base_idx); cudaFree(ws->d_spread_weight);
    cudaFree(ws->d_shift_r); cudaFree(ws->d_shift_i);
    cudaFree(ws->d_fft_real); cudaFree(ws->d_fft_imag);
    cudaFree(ws->d_y_real); cudaFree(ws->d_y_imag);
    cudaFree(ws->d_out_real); cudaFree(ws->d_out_imag);
    free(ws);
}

void cunufft_free_plan(cunufft_plan *p) {
    if (!p) return;
    cudaFree(p->d_deconv);
    cudaFree(p->d_twiddle_real); cudaFree(p->d_twiddle_imag);
    free(p);
}
