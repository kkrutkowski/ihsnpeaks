/*
 * lc_model.c — Phased model overlay: single-frequency model evaluation.
 *
 * Split from lc_periodogram.c for maintainability. This file is #include'd at the
 * end of lc_periodogram.c (single-TU compilation) because the upstream header-only
 * libraries (sds.h, convolution.h) contain non-static definitions that must only be
 * compiled once.
 *
 * Computes model y-values at each data point for the pivot frequency.
 * Method-specific:
 *   IHS  -> direct sine/cosine sums (skip Szego orthogonalization)
 *   AoV  -> Szego recursion trig polynomial
 *   GB   -> convolution smoother (scatter style)
 *   BLS  -> boxcar fit (line style)
 *
 * AoV and IHS use GNU vector extensions with custom vectorized sin/cos
 * (from trig.h / nufft1.c) and zero-pad buffers to the vector length so
 * all inner loops are fully SIMD with no scalar remainders.
 */

/* ===================================================================================
 * Vectorized trig helpers (GNU vector extensions)
 * =================================================================================== */

#define LC_VECF_LEN (VEC_BYTES / 4)
typedef float lc_vecf __attribute__((vector_size(VEC_BYTES)));
typedef int32_t lc_veci __attribute__((vector_size(VEC_BYTES)));

/* Element-wise vector select: where mask is true (all-ones), pick a; else pick b.
 * GCC vector extensions do not support the ternary operator on vector types,
 * so we use bitwise AND/OR with the comparison mask. */
static inline lc_vecf lc_vselect(lc_veci mask, lc_vecf a, lc_vecf b) {
    return (lc_vecf)(((lc_veci)a & mask) | ((lc_veci)b & ~mask));
}

/* Vectorized sin(2*pi*x) — replicates sin2pif_tls polynomial approximation
 * element-wise on a GNU vector type. Range reduction uses vector comparisons
 * compiled to blends/masks by the backend. */
static inline lc_vecf lc_vecf_sin2pi(lc_vecf x) {
    lc_vecf zero = (lc_vecf){};
    lc_vecf one = zero + 1.0f;
    lc_vecf half = zero + 0.5f;
    lc_vecf quarter = zero + 0.25f;

    /* Fractional part */
    lc_veci xi = __builtin_convertvector(x, lc_veci);
    lc_vecf fi = __builtin_convertvector(xi, lc_vecf);
    lc_vecf f = x - fi;
    f = lc_vselect(f < zero, f + one, f);

    /* Sign and fold */
    lc_vecf sign = one;
    lc_veci ge_half = f >= half;
    sign = lc_vselect(ge_half, zero - one, sign);
    f = lc_vselect(ge_half, f - half, f);
    lc_veci gt_quarter = f > quarter;
    f = lc_vselect(gt_quarter, half - f, f);

    /* Polynomial (Horner) */
    lc_vecf f2 = f * f;
    lc_vecf p = f2 * (zero + 39.536706065730207835108712734262f) + (zero - 76.549782293595742666226937116116f);
    p = p * f2 + (zero + 81.601004073261773523492199897936f);
    p = p * f2 + (zero - 41.341655031416278077153126232486f);
    p = p * f2 + (zero + 6.2831851600894774430188071795666f);
    p = p * f;
    return p * sign;
}

/* Vectorized cos(2*pi*x) — replicates cos2pif_tls polynomial approximation. */
static inline lc_vecf lc_vecf_cos2pi(lc_vecf x) {
    lc_vecf zero = (lc_vecf){};
    lc_vecf one = zero + 1.0f;
    lc_vecf half = zero + 0.5f;
    lc_vecf quarter = zero + 0.25f;

    lc_veci xi = __builtin_convertvector(x, lc_veci);
    lc_vecf fi = __builtin_convertvector(xi, lc_vecf);
    lc_vecf f = x - fi;
    f = lc_vselect(f < zero, f + one, f);
    lc_veci gt_half = f > half;
    f = lc_vselect(gt_half, one - f, f);

    lc_vecf sign = one;
    lc_veci gt_quarter = f > quarter;
    sign = lc_vselect(gt_quarter, zero - one, sign);
    f = lc_vselect(gt_quarter, half - f, f);

    lc_vecf f2 = f * f;
    lc_vecf p = f2 * (zero + 56.242380464873243259663276802701f) + (zero - 85.240330322699427859509454517828f);
    p = p * f2 + (zero + 64.934590626780991246193352727536f);
    p = p * f2 + (zero - 19.739171434702393618770795066531f);
    p = p * f2 + (zero + 0.99999995346667013630639784578184f);
    return p * sign;
}

/* Round up to the next multiple of LC_VECF_LEN (for zero-padding). */
static inline uint32_t lc_model_pad_len(uint32_t n) {
    return (n + LC_VECF_LEN - 1U) & ~((uint32_t)LC_VECF_LEN - 1U);
}

/* -----------------------------------------------------------------------------------
 * AoV model: Szego recursion at a single frequency (vectorized inner loop)
 * Reference: aovdist (aovmh routine)
 *
 * The reference multiplies data by exp(i*nh*angle) before projecting onto Szego
 * polynomials (cfr/cfi initialization). The model reconstruction must undo this
 * frequency shift by multiplying back by exp(-i*nh*angle):
 *   model[k] += Re(coeff * P_n(z[k]) * exp(-i*nh*angle[k]))
 *            = A*cos(nh*angle[k]) + B*sin(nh*angle[k])
 * where A = Re(coeff*P_n), B = Im(coeff*P_n).
 * ----------------------------------------------------------------------------------- */
static int lc_model_aov(const double *x, const float *y, const float *dy, uint32_t n, double freq, int nterms, float *model_out) {
    if (n == 0 || nterms < 1) return -1;
    int nh = nterms;
    int nn2 = nh + nh;

    if ((int)n < nn2 + 2) return -1; /* not enough DOF */

    uint32_t npad = lc_model_pad_len(n);

    /* Allocate zero-padded arrays (aligned for SIMD loads/stores) */
    float *zr = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *zi = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *pr = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *pi = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *znr = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *zni = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *cfr = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *cfi = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *v = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *rw = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *model_f = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *phase_nh = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    if (!zr || !zi || !pr || !pi || !znr || !zni || !cfr || !cfi || !v || !rw || !model_f || !phase_nh) {
        free(zr); free(zi); free(pr); free(pi); free(znr); free(zni);
        free(cfr); free(cfi); free(v); free(rw); free(model_f); free(phase_nh);
        return -1;
    }
    memset(zr, 0, (size_t)npad * sizeof(float));
    memset(zi, 0, (size_t)npad * sizeof(float));
    memset(pr, 0, (size_t)npad * sizeof(float));
    memset(pi, 0, (size_t)npad * sizeof(float));
    memset(znr, 0, (size_t)npad * sizeof(float));
    memset(zni, 0, (size_t)npad * sizeof(float));
    memset(cfr, 0, (size_t)npad * sizeof(float));
    memset(cfi, 0, (size_t)npad * sizeof(float));
    memset(v, 0, (size_t)npad * sizeof(float));
    memset(rw, 0, (size_t)npad * sizeof(float));
    memset(model_f, 0, (size_t)npad * sizeof(float));
    memset(phase_nh, 0, (size_t)npad * sizeof(float));

    /* Normalize: subtract weighted mean (or unweighted mean if dy is NULL) */
    float sav = 0.0f;
    float ws = 0.0f;
    float epsilon = 1e-3f; /* default parameters epsilon */
    if (dy) {
        for (uint32_t k = 0; k < n; k++) {
            float weight = 1.0f / (dy[k] * dy[k] + epsilon);
            rw[k] = sqrtf(weight);
            sav += y[k] * weight;
            ws += weight;
        }
        sav /= ws;
        for (uint32_t k = 0; k < n; k++) {
            v[k] = (y[k] - sav) * rw[k];
        }
    } else {
        for (uint32_t k = 0; k < n; k++) {
            sav += y[k];
            rw[k] = 1.0f;
        }
        sav /= (float)n;
        for (uint32_t k = 0; k < n; k++) {
            v[k] = y[k] - sav;
        }
        for (uint32_t k = n; k < npad; k++) {
            rw[k] = 1.0f;
        }
    }

    /* Compute phases in DOUBLE precision, then initialize Szego arrays.
     * Phase = fmod((x[k] - x[0]) * freq, 1.0) — double precision as required. */
    double t0 = x[0];
    for (uint32_t k = 0; k < n; k++) {
        double dph = (x[k] - t0) * freq;
        dph -= floor(dph);
        float phase_f = (float)dph; /* single precision for trig fitting */
        zr[k] = cos2pif_tls(phase_f);
        zi[k] = sin2pif_tls(phase_f);
        pr[k] = rw[k];
        pi[k] = 0.0f;
        znr[k] = 1.0f;
        zni[k] = 0.0f;
        float angle_nh = phase_f * (float)nh;
        phase_nh[k] = angle_nh; /* store for de-rotation during model accumulation */
        cfr[k] = v[k] * cos2pif_tls(angle_nh);
        cfi[k] = v[k] * sin2pif_tls(angle_nh);
    }
    for (uint32_t k = n; k < npad; k++) {
        znr[k] = 1.0f;
    }

    /* Szego recursion with SIMD-vectorized inner loop.
     * Accumulate model = sum_n Re(coeff_n * P_n * exp(-i*nh*angle)) at each step. */
    for (int step = 0; step <= nn2; step++) {
        /* Reductions: sn, alr, ali, scr, sci — vectorized accumulation */
        lc_vecf v_sn = (lc_vecf){};
        lc_vecf v_alr = (lc_vecf){};
        lc_vecf v_ali = (lc_vecf){};
        lc_vecf v_scr = (lc_vecf){};
        lc_vecf v_sci = (lc_vecf){};

        for (uint32_t k = 0; k < npad; k += LC_VECF_LEN) {
            lc_vecf vpr = *(lc_vecf *)&pr[k];
            lc_vecf vpi = *(lc_vecf *)&pi[k];
            lc_vecf vzr = *(lc_vecf *)&zr[k];
            lc_vecf vzi = *(lc_vecf *)&zi[k];
            lc_vecf vcfr = *(lc_vecf *)&cfr[k];
            lc_vecf vcfi = *(lc_vecf *)&cfi[k];
            lc_vecf vrw = *(lc_vecf *)&rw[k];

            v_sn += vpr * vpr + vpi * vpi;
            v_alr += (vzr * vpr - vzi * vpi) * vrw;
            v_ali += (vzr * vpi + vzi * vpr) * vrw;
            v_scr += vpr * vcfr + vpi * vcfi;
            v_sci += vpr * vcfi - vpi * vcfr;
        }

        /* Horizontal sum */
        float sn = 0.0f, alr = 0.0f, ali = 0.0f, scr = 0.0f, sci = 0.0f;
        for (int lane = 0; lane < LC_VECF_LEN; lane++) {
            sn += v_sn[lane];
            alr += v_alr[lane];
            ali += v_ali[lane];
            scr += v_scr[lane];
            sci += v_sci[lane];
        }

        if (sn < 1e-30f) sn = 1e-30f;
        alr /= sn;
        ali /= sn;

        /* Fourier coefficient at this Szego step */
        float coeff_r = scr / sn;
        float coeff_i = sci / sn;

        /* Accumulate model with de-rotation by exp(-i*nh*angle):
         * A = coeff_r*pr - coeff_i*pi  (Re of coeff * P_n)
         * B = coeff_r*pi + coeff_i*pr  (Im of coeff * P_n)
         * model[k] += A*cos(nh*angle[k]) + B*sin(nh*angle[k]) */
        lc_vecf vcr = (lc_vecf){} + coeff_r;
        lc_vecf vci = (lc_vecf){} + coeff_i;
        for (uint32_t k = 0; k < npad; k += LC_VECF_LEN) {
            lc_vecf vm = *(lc_vecf *)&model_f[k];
            lc_vecf vpr = *(lc_vecf *)&pr[k];
            lc_vecf vpi = *(lc_vecf *)&pi[k];
            lc_vecf vph = *(lc_vecf *)&phase_nh[k];
            lc_vecf vc = lc_vecf_cos2pi(vph);
            lc_vecf vs = lc_vecf_sin2pi(vph);
            lc_vecf vA = vcr * vpr - vci * vpi;
            lc_vecf vB = vcr * vpi + vci * vpr;
            *(lc_vecf *)&model_f[k] = vm + vA * vc + vB * vs;
        }

        /* Szego polynomial update (vectorized) */
        lc_vecf valr = (lc_vecf){} + alr;
        lc_vecf vali = (lc_vecf){} + ali;
        for (uint32_t k = 0; k < npad; k += LC_VECF_LEN) {
            lc_vecf vpr = *(lc_vecf *)&pr[k];
            lc_vecf vpi = *(lc_vecf *)&pi[k];
            lc_vecf vzr = *(lc_vecf *)&zr[k];
            lc_vecf vzi = *(lc_vecf *)&zi[k];
            lc_vecf vznr = *(lc_vecf *)&znr[k];
            lc_vecf vzni = *(lc_vecf *)&zni[k];
            lc_vecf vrw = *(lc_vecf *)&rw[k];

            /* sr = alr*znr - ali*zni; si = alr*zni + ali*znr */
            lc_vecf vsr = valr * vznr - vali * vzni;
            lc_vecf vsi = valr * vzni + vali * vznr;

            /* new pr = pr*zr - pi*zi - sr*pr - si*pi */
            lc_vecf new_pr = vpr * vzr - vpi * vzi - vsr * vpr - vsi * vpi;
            /* new pi = pr*zi + pi*zr + sr*pi - si*pr */
            lc_vecf new_pi = vpr * vzi + vpi * vzr + vsr * vpi - vsi * vpr;

            *(lc_vecf *)&pr[k] = new_pr;
            *(lc_vecf *)&pi[k] = new_pi;

            /* znr/zni advance: multiply by z */
            lc_vecf new_znr = vznr * vzr - vzni * vzi;
            lc_vecf new_zni = vzni * vzr + vznr * vzi;
            *(lc_vecf *)&znr[k] = new_znr;
            *(lc_vecf *)&zni[k] = new_zni;
        }
    }

    /* Copy model (add mean back, divide by weight factor if weighted) */
    for (uint32_t k = 0; k < n; k++) {
        float val = model_f[k];
        if (dy && rw[k] > 0.0f) {
            val /= rw[k];
        }
        model_out[k] = val + sav;
    }

    free(zr); free(zi); free(pr); free(pi); free(znr); free(zni);
    free(cfr); free(cfi); free(v); free(rw); free(model_f); free(phase_nh);
    return 0;
}

/* -----------------------------------------------------------------------------------
 * IHS model: direct sine/cosine sums (skip Szego orthogonalization).
 * Assumes measurements are already orthogonal — uses simple trig projection.
 * SIMD-vectorized with zero-padded buffers.
 * ----------------------------------------------------------------------------------- */
static int lc_model_ihs(const double *x, const float *y, uint32_t n, double freq, int nterms, float *model_out) {
    if (n == 0 || nterms < 1) return -1;

    uint32_t npad = lc_model_pad_len(n);

    float *phase_f = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *v = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    float *model_f = (float *)aligned_alloc(64, (size_t)npad * sizeof(float));
    if (!phase_f || !v || !model_f) {
        free(phase_f); free(v); free(model_f);
        return -1;
    }
    memset(phase_f, 0, (size_t)npad * sizeof(float));
    memset(v, 0, (size_t)npad * sizeof(float));
    memset(model_f, 0, (size_t)npad * sizeof(float));

    /* Mean-subtract */
    float mean = 0.0f;
    for (uint32_t k = 0; k < n; k++) mean += y[k];
    mean /= (float)n;
    for (uint32_t k = 0; k < n; k++) v[k] = y[k] - mean;

    /* Compute phases in DOUBLE precision */
    double t0 = x[0];
    for (uint32_t k = 0; k < n; k++) {
        double dph = (x[k] - t0) * freq;
        dph -= floor(dph);
        phase_f[k] = (float)dph;
    }

    /* For each harmonic q, compute Cq and Sq via vectorized dot products,
     * then accumulate model += Cq*cos(2*pi*q*phase) + Sq*sin(2*pi*q*phase). */
    float norm = 2.0f / (float)n;

    for (int q = 1; q <= nterms; q++) {
        float qf = (float)q;

        /* Vectorized accumulation of Cq = sum(v*cos), Sq = sum(v*sin) */
        lc_vecf v_cq = (lc_vecf){};
        lc_vecf v_sq = (lc_vecf){};

        for (uint32_t k = 0; k < npad; k += LC_VECF_LEN) {
            lc_vecf vp = *(lc_vecf *)&phase_f[k];
            lc_vecf angle = vp * ((lc_vecf){} + qf);
            lc_vecf vc = lc_vecf_cos2pi(angle);
            lc_vecf vs = lc_vecf_sin2pi(angle);
            lc_vecf vv = *(lc_vecf *)&v[k];
            v_cq += vv * vc;
            v_sq += vv * vs;
        }

        float Cq = 0.0f, Sq = 0.0f;
        for (int lane = 0; lane < LC_VECF_LEN; lane++) {
            Cq += v_cq[lane];
            Sq += v_sq[lane];
        }
        Cq *= norm;
        Sq *= norm;

        /* Accumulate model: model[k] += Cq*cos(2*pi*q*phase[k]) + Sq*sin(2*pi*q*phase[k]) */
        lc_vecf vCq = (lc_vecf){} + Cq;
        lc_vecf vSq = (lc_vecf){} + Sq;
        for (uint32_t k = 0; k < npad; k += LC_VECF_LEN) {
            lc_vecf vp = *(lc_vecf *)&phase_f[k];
            lc_vecf angle = vp * ((lc_vecf){} + qf);
            lc_vecf vc = lc_vecf_cos2pi(angle);
            lc_vecf vs = lc_vecf_sin2pi(angle);
            lc_vecf vm = *(lc_vecf *)&model_f[k];
            *(lc_vecf *)&model_f[k] = vm + vCq * vc + vSq * vs;
        }
    }

    /* Copy model (add mean back) */
    for (uint32_t k = 0; k < n; k++) {
        model_out[k] = model_f[k] + mean;
    }

    free(phase_f); free(v); free(model_f);
    return 0;
}

/* -----------------------------------------------------------------------------------
 * GB model: convolution smoother at a single frequency.
 * The convolution-based model is not continuous by default, so it is rendered
 * as a scatter plot.
 * TODO: FFT interpolation via zero-padded FFT/IFFT + shift theorem could produce
 * a continuous model curve. Not implemented yet.
 * ----------------------------------------------------------------------------------- */
static int lc_model_gb(const double *x, const float *y, const float *dy, uint32_t n, double freq, float gbAlpha, float *model_out) {
    if (n == 0) return -1;

    /* Allocate a minimal buffer_t for gb_prepare_projection */
    buffer_t buf;
    memset(&buf, 0, sizeof(buf));
    buf.n = n;
    buf.x = (double *)x;
    buf.y = (float *)y;
    buf.dy = (float *)dy;

    /* Allocate sort/convolution scratch */
    size_t kv_bytes = (size_t)n * sizeof(kvpair);
    size_t idx_bytes = 1024 * sizeof(size_t);
    void *buf0 = calloc(1, kv_bytes > idx_bytes ? kv_bytes : idx_bytes);
    void *buf1 = calloc(1, kv_bytes);
    void *buf2 = calloc(1, (size_t)n * sizeof(double));
    size_t *pidx = (size_t *)calloc(1024, sizeof(size_t));
    if (!buf0 || !buf1 || !buf2 || !pidx) {
        free(buf0); free(buf1); free(buf2); free(pidx);
        return -1;
    }

    void *bufs[3] = {buf0, buf1, buf2};
    buf.buf = bufs;
    buf.pidx = pidx;

    gb_projection_t projection = {0};
    if (!gb_prepare_projection(&buf, freq, true, gbAlpha, 1.0, &projection)) {
        free(buf0); free(buf1); free(buf2); free(pidx);
        return -1;
    }

    /* Map model values back to original indices.
     * model[i] in projection is the smoothed value for sorted[i];
     * the fitted model = scale * smoothed. */
    double scale = clamp_gbls_scale(projection.scale);
    for (int i = 0; i < projection.n; i++) {
        uint32_t orig_idx = projection.sorted[i].parts.idx;
        if (orig_idx < n) {
            model_out[orig_idx] = (float)(scale * projection.model[i]);
        }
    }

    free(buf0); free(buf1); free(buf2); free(pidx);
    return 0;
}

/* -----------------------------------------------------------------------------------
 * BLS model: boxcar fit at a single frequency.
 * Passes the same options (blsMinRelWidth, blsMaxRelWidth, blsWidthCount) as
 * the period search and returns the fitted boxcar levels.
 * ----------------------------------------------------------------------------------- */
static int lc_model_bls(const double *x, const float *y, const float *dy, uint32_t n, double freq,
                        double blsMinRelWidth, double blsMaxRelWidth, int blsWidthCount, float *model_out) {
    if (n <= 2) return -1;

    /* Allocate a minimal buffer_t for get_bls_result */
    buffer_t buf;
    memset(&buf, 0, sizeof(buf));
    buf.n = n;
    buf.x = (double *)x;
    buf.y = (float *)y;
    buf.dy = (float *)dy;

    size_t idx_bytes = 1024 * sizeof(size_t);
    size_t sort_bytes = (size_t)n * sizeof(uint64_t);
    void *buf0 = calloc(1, sort_bytes > idx_bytes ? sort_bytes : idx_bytes);
    void *buf1 = calloc(1, sort_bytes);
    size_t *pidx = (size_t *)calloc(1024, sizeof(size_t));
    if (!buf0 || !buf1 || !pidx) {
        free(buf0); free(buf1); free(pidx);
        return -1;
    }

    void *bufs[2] = {buf0, buf1};
    buf.buf = bufs;
    buf.pidx = pidx;

    /* Build parameters for BLS */
    parameters params;
    memset(&params, 0, sizeof(params));
    params.blsMinRelWidth = blsMinRelWidth > 0.0 ? blsMinRelWidth : 0.01;
    params.blsMaxRelWidth = blsMaxRelWidth > 0.0 ? blsMaxRelWidth : 0.5;
    params.blsWidthCount = blsWidthCount > 0 ? blsWidthCount : 10;

    /* Replicate get_bls_result logic to extract in/out levels and positions. */
    int nn = (int)n;
    uint64_t *input = (uint64_t *)buf0;
    uint64_t *sorted = (uint64_t *)buf1;
    uint64_t **sort_in = (uint64_t **)(&bufs[0]);
    uint64_t **sort_out = (uint64_t **)(&bufs[1]);

    bls_sums_t total = {0};
    for (int i = 0; i < nn; ++i) {
        uint16_t key = phase_key_10(x[i] * freq);
        input[i] = bls_pack_phase_index(key, (uint32_t)i);
        bls_add_index(&buf, (uint32_t)i, 1.0, &total);
    }
    csort64_10(sort_in, buf.n, sort_out, buf.pidx);

    double ref_sse = total.y2 - ((total.y * total.y) / (double)nn);
    if (!(ref_sse > 0.0)) {
        free(buf0); free(buf1); free(pidx);
        return -1;
    }

    int width_count = params.blsWidthCount;
    double best_score = -1.0;
    double best_in_level = 0.0;
    double best_out_level = 0.0;
    int best_start = -1;
    int best_width = 0;

    for (int width_idx = 0; width_idx < width_count; ++width_idx) {
        double rel_width = bls_relative_width_at(&params, width_idx);
        if (!(rel_width > 0.0)) continue;
        int width = (int)ceil(rel_width * (double)nn);
        if (width < 1) width = 1;
        if (width >= nn) width = nn - 1;

        bls_sums_t in = {0};
        for (int j = 0; j < width; ++j) {
            bls_add_sorted(&buf, sorted, j, 1.0, &in);
        }

        int out_count = nn - width;
        for (int start = 0; start < nn; ++start) {
            bls_sums_t out = {.y = total.y - in.y, .y2 = total.y2 - in.y2};
            if (out_count > 0) {
                double in_sse = in.y2 - ((in.y * in.y) / (double)width);
                double out_sse = out.y2 - ((out.y * out.y) / (double)out_count);
                if (in_sse < 0.0) in_sse = 0.0;
                if (out_sse < 0.0) out_sse = 0.0;
                double model_sse = in_sse + out_sse;
                double explained = ref_sse - model_sse;
                if (explained < 0.0) explained = 0.0;
                if (explained > best_score) {
                    best_score = explained;
                    best_in_level = in.y / (double)width;
                    best_out_level = out.y / (double)out_count;
                    best_start = start;
                    best_width = width;
                }
            }
            if (start + 1 < nn) {
                bls_add_sorted(&buf, sorted, start, -1.0, &in);
                bls_add_sorted(&buf, sorted, start + width, 1.0, &in);
            }
        }
    }

    if (best_start < 0) {
        free(buf0); free(buf1); free(pidx);
        return -1;
    }

    /* Mark which original indices fall in the transit window */
    uint8_t *in_transit = (uint8_t *)calloc(n, 1);
    if (!in_transit) {
        free(buf0); free(buf1); free(pidx);
        return -1;
    }
    for (int j = 0; j < best_width; ++j) {
        int pos = best_start + j;
        if (pos >= nn) pos %= nn;
        uint32_t orig_idx = bls_unpack_index(sorted[pos]);
        if (orig_idx < n) in_transit[orig_idx] = 1;
    }

    for (uint32_t i = 0; i < n; i++) {
        model_out[i] = (float)(in_transit[i] ? best_in_level : best_out_level);
    }

    free(in_transit);
    free(buf0); free(buf1); free(pidx);
    return 0;
}

/* -----------------------------------------------------------------------------------
 * Main dispatch: lc_compute_phased_model
 * ----------------------------------------------------------------------------------- */
LC_API int lc_compute_phased_model(const lc_data_t *data, const lc_periodogram_config_t *cfg, double freq, float *model_out,
                                   lc_model_style_t *style_out) {
    if (!data || data->n == 0 || !cfg || !model_out || !style_out) return -1;
    if (!(freq > 0.0)) return -1;

    uint32_t n = data->n;
    int nterms = cfg->nterms > 0 ? cfg->nterms : 3;

    /* Set default style */
    *style_out = LC_MODEL_LINE;

    switch (cfg->method) {
    case LC_SPEC_IHS:
        return lc_model_ihs(data->x, data->y, n, freq, nterms, model_out);

    case LC_SPEC_AOV:
        return lc_model_aov(data->x, data->y, data->dy, n, freq, nterms, model_out);

    case LC_SPEC_GB: {
        /* GB: convolution smoother. Scatter style (not continuous).
         * The period search works on zero-mean preprocessed data, so we subtract
         * the mean before fitting and add it back afterwards. */
        *style_out = LC_MODEL_SCATTER;
        float gbAlpha = (float)(cfg->oversmoothing > 0.0 ? cfg->oversmoothing : 0.2);
        float mean = 0.0f;
        for (uint32_t i = 0; i < n; i++) mean += data->y[i];
        mean /= (float)n;
        float *y_centered = (float *)malloc((size_t)n * sizeof(float));
        if (!y_centered) return -1;
        for (uint32_t i = 0; i < n; i++) y_centered[i] = data->y[i] - mean;
        int rc = lc_model_gb(data->x, y_centered, data->dy, n, freq, gbAlpha, model_out);
        free(y_centered);
        if (rc != 0) return rc;
        for (uint32_t i = 0; i < n; i++) model_out[i] += mean;
        return 0;
    }

    case LC_SPEC_BLS: {
        /* BLS: boxcar fit. Line style.
         * Same mean-handling as GB: subtract before fit, add back after. */
        double blsMinRelWidth = cfg->oversmoothing > 0.0 ? cfg->oversmoothing : 0.01;
        double blsMaxRelWidth = 0.5;
        int blsWidthCount = cfg->nbins > 0 ? cfg->nbins : 10;
        float mean = 0.0f;
        for (uint32_t i = 0; i < n; i++) mean += data->y[i];
        mean /= (float)n;
        float *y_centered = (float *)malloc((size_t)n * sizeof(float));
        if (!y_centered) return -1;
        for (uint32_t i = 0; i < n; i++) y_centered[i] = data->y[i] - mean;
        int rc = lc_model_bls(data->x, y_centered, data->dy, n, freq, blsMinRelWidth, blsMaxRelWidth, blsWidthCount, model_out);
        free(y_centered);
        if (rc != 0) return rc;
        for (uint32_t i = 0; i < n; i++) model_out[i] += mean;
        return 0;
    }

    default:
        return -1;
    }
}
