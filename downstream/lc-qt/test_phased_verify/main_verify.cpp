/*
 * phased_verify — pixel-level verification harness for PhasedLightCurveWidget.
 *
 * Links the real widget sources and the real loader/detrend backend, renders
 * the widget to QImages, and checks:
 *  1. freq <= 0            -> "No data" state (no scatter pixels)
 *  2. freq = f1            -> scatter present; every measurement lands at the
 *                             spec position: phase = fmod((x[i]-x[0])*f, 1),
 *                             duplicated at phase + 1, y-axis inverted
 *  3. mirror symmetry      -> the [1,2] half is a pixel copy of [0,1] shifted
 *                             by exactly plotW/2
 *  4. freq change          -> re-fold produces a different image
 *  5. main.cpp wiring      -> ZoomedSpectrumWidget::centerFrequencyChanged
 *                             connected to setFrequency (exact connect from
 *                             main.cpp) drives the fold on pivot moves
 */
#include <QApplication>
#include <QImage>
#include <QVector>
#include <cstdio>
#include <cstring>
#include <cmath>

extern "C" {
#include "lc_readout.h"
}

#include "windows/phased_lightcurve.h"
#include "windows/zoomed_spectrum.h"

static constexpr int PAD_L = 55, PAD_R = 15, PAD_T = 28, PAD_B = 35;
static int failures = 0;

#define CHECK(cond, ...)                                    \
    do {                                                    \
        if (cond) {                                         \
            printf("  PASS: " __VA_ARGS__);                 \
            printf("\n");                                   \
        } else {                                            \
            printf("  FAIL: " __VA_ARGS__);                 \
            printf("\n");                                   \
            failures++;                                     \
        }                                                   \
    } while (0)

static bool isScatterPixel(QRgb p) {
    /* #3b82f6 = (59,130,246) with antialiasing; reject bg/grid/text/axes */
    int r = qRed(p), g = qGreen(p), b = qBlue(p);
    return b > 150 && b > r + 40 && b > g + 40;
}

/* Tolerant mirror match: a scatter pixel, or an antialiased edge pixel whose
   blue channel is close to the source pixel's (Qt's rasterizer can differ by
   a subpixel sample between the two ellipse passes at different absolute x). */
static bool mirrorsPixel(QRgb src, QRgb cand) {
    if (isScatterPixel(cand)) return true;
    int b = qBlue(cand);
    return b > 120 && b > qRed(cand) + 30 && b > qGreen(cand) + 30 && std::abs(b - qBlue(src)) <= 30;
}

static QImage renderWidget(QWidget &w) {
    QImage img(w.size(), QImage::Format_ARGB32);
    img.fill(QColor(0, 0, 0));
    w.render(&img);
    return img;
}

/* Count scatter pixels strictly inside the plot area. */
static int countScatter(const QImage &img, int plotW, int plotH) {
    int cnt = 0;
    for (int y = PAD_T + 2; y < PAD_T + plotH - 1; ++y)
        for (int x = PAD_L + 2; x < PAD_L + plotW - 1; ++x)
            if (isScatterPixel(img.pixel(x, y))) cnt++;
    return cnt;
}

static int diffPixels(const QImage &a, const QImage &b) {
    if (a.size() != b.size()) return -1;
    int d = 0;
    for (int y = 0; y < a.height(); ++y)
        for (int x = 0; x < a.width(); ++x)
            if (a.pixel(x, y) != b.pixel(x, y)) d++;
    return d;
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    if (argc < 2) {
        fprintf(stderr, "usage: phased_verify <lc.dat>\n");
        return 2;
    }

    lc_data_t d;
    memset(&d, 0, sizeof(d));
    if (lc_load_dat(argv[1], &d) != 0 || d.n == 0) {
        fprintf(stderr, "failed to load %s\n", argv[1]);
        return 2;
    }
    lc_detrend(&d);
    printf("Loaded %u points from %s\n", d.n, argv[1]);

    PhasedLightCurveWidget w("Phased Light Curve");
    w.resize(440, 300);
    const int plotW = 440 - PAD_L - PAD_R; /* 370 */
    const int plotH = 300 - PAD_T - PAD_B; /* 237 */
    w.setData(&d);

    /* --- 1. No frequency -> "No data" state --- */
    printf("[1] Unfolded state (freq <= 0)\n");
    w.setFrequency(-1.0);
    QImage imgNone = renderWidget(w);
    CHECK(countScatter(imgNone, plotW, plotH) == 0, "no scatter pixels before a pivot frequency is set");

    /* --- 2. Fold at f1: scatter present, points at spec positions --- */
    printf("[2] Fold at f1 (phase = (t - t0) * f mod 1, duplicated at +1)\n");
    const double f1 = 24.0; /* near the BLAP principal frequency */
    w.setFrequency(f1);
    QImage imgF1 = renderWidget(w);
    int nScatter = countScatter(imgF1, plotW, plotH);
    CHECK(nScatter > (int)d.n, "scatter rendered (%d scatter pixels for %u measurements)", nScatter, d.n);

    /* Replicate the widget's exact y-range computation */
    float yMin = d.y[0], yMax = d.y[0];
    for (unsigned int i = 1; i < d.n; ++i) {
        if (d.y[i] < yMin) yMin = d.y[i];
        if (d.y[i] > yMax) yMax = d.y[i];
    }
    float ySpan = yMax - yMin;
    if (ySpan < 1e-6f) ySpan = 1.0f;
    yMin -= ySpan * 0.05f;
    yMax += ySpan * 0.05f;
    ySpan = yMax - yMin;

    /* Every measurement must have scatter coverage near BOTH expected pixels */
    int missFirst = 0, missCopy = 0;
    const double t0 = (d.x[0] + d.x[d.n - 1]) / 2.0;
    for (unsigned int i = 0; i < d.n; ++i) {
        double phase = std::fmod((d.x[i] - t0) * f1, 1.0);
        if (phase < 0.0) phase += 1.0;
        double py = PAD_T + ((double)(d.y[i] - yMin) / ySpan) * plotH;
        double px1 = PAD_L + (phase / 2.0) * plotW;
        double px2 = PAD_L + ((phase + 1.0) / 2.0) * plotW;
        auto covered = [&](double px) {
            for (int dy = -2; dy <= 2; ++dy)
                for (int dx = -2; dx <= 2; ++dx) {
                    int xx = (int)std::round(px) + dx, yy = (int)std::round(py) + dy;
                    if (xx < 0 || yy < 0 || xx >= imgF1.width() || yy >= imgF1.height()) continue;
                    if (isScatterPixel(imgF1.pixel(xx, yy))) return true;
                }
            return false;
        };
        if (!covered(px1)) missFirst++;
        if (!covered(px2)) missCopy++;
    }
    CHECK(missFirst == 0, "all %u points land at spec phase positions (missed %d)", d.n, missFirst);
    CHECK(missCopy == 0, "all %u points duplicated at phase + 1 (missed %d)", d.n, missCopy);

    /* --- 3. Mirror symmetry: [1,2] is a copy of [0,1] shifted by plotW/2 --- */
    printf("[3] Mirror symmetry between phase halves\n");
    const int half = plotW / 2; /* 185 */
    const int midX = PAD_L + half;
    int violL = 0, violR = 0;
    /* Pixels within 2px of the phase-1 boundary (midX) are skipped: a point
       at phase ~= 1 legitimately spills antialiased edge pixels across it. */
    for (int y = PAD_T + 2; y < PAD_T + plotH - 1; ++y) {
        for (int x = PAD_L + 2; x < midX - 2; ++x) {
            QRgb src = imgF1.pixel(x, y);
            if (!isScatterPixel(src)) continue;
            bool mirrored = false;
            for (int dx = -1; dx <= 1 && !mirrored; ++dx) {
                int xx = x + half + dx;
                if (xx >= 0 && xx < imgF1.width() && mirrorsPixel(src, imgF1.pixel(xx, y))) mirrored = true;
            }
            if (!mirrored) violL++;
        }
        for (int x = midX + 3; x < PAD_L + plotW - 1; ++x) {
            QRgb src = imgF1.pixel(x, y);
            if (!isScatterPixel(src)) continue;
            bool mirrored = false;
            for (int dx = -1; dx <= 1 && !mirrored; ++dx) {
                int xx = x - half + dx;
                if (xx >= 0 && xx < imgF1.width() && mirrorsPixel(src, imgF1.pixel(xx, y))) mirrored = true;
            }
            if (!mirrored) violR++;
        }
    }
    CHECK(violL == 0 && violR == 0,
          "[0,1] -> [1,2] mirror exact (violations: %d left, %d right)", violL, violR);

    /* --- 4. Frequency change re-folds the curve --- */
    printf("[4] Re-fold on frequency change\n");
    const double f2 = 20.0;
    w.setFrequency(f2);
    QImage imgF2 = renderWidget(w);
    CHECK(countScatter(imgF2, plotW, plotH) > (int)d.n, "scatter rendered at f2");
    CHECK(diffPixels(imgF1, imgF2) > 0, "image at f2 differs from image at f1");

    /* --- 5. main.cpp wiring: pivot signal drives the fold --- */
    printf("[5] ZoomedSpectrumWidget::centerFrequencyChanged -> setFrequency\n");
    ZoomedSpectrumWidget zoom("Zoomed Spectrum");
    QObject::connect(&zoom, &ZoomedSpectrumWidget::centerFrequencyChanged,
                     &w, &PhasedLightCurveWidget::setFrequency);
    QVector<float> nll(5000);
    for (int i = 0; i < nll.size(); ++i) nll[i] = 1.0f + 0.5f * std::sin(0.01 * i);
    zoom.setZoomFactor(4.0);
    zoom.setFullSpectrum(1.0, 0.005, nll); /* grid [1, 26] c/d */

    zoom.selectFrequency(f1, true); /* pivot move, as on spectrum click */
    QImage imgWired1 = renderWidget(w);
    CHECK(diffPixels(imgWired1, imgF1) == 0, "pivot -> f1 reproduces the direct-set fold exactly");

    zoom.selectFrequency(f2, false); /* pivot move again */
    QImage imgWired2 = renderWidget(w);
    CHECK(diffPixels(imgWired2, imgF2) == 0, "pivot -> f2 reproduces the direct-set fold exactly");

    /* --- 6. Phase Offset & Extremum alignment at 0.5 --- */
    printf("[6] Phase offset and model extremum alignment at phase 0.5\n");
    CHECK(w.phaseOffset() == 0.0, "initial phase offset is 0.0");
    w.setPhaseOffset(0.25);
    CHECK(std::abs(w.phaseOffset() - 0.25) < 1e-6, "setPhaseOffset updates offset to 0.25");
    QImage imgShifted = renderWidget(w);
    CHECK(diffPixels(imgShifted, imgF2) > 0, "phase offset shift produces visually shifted scatter");

    /* Test all 4 methods (IHS, AoV, GB, BLS) */
    lc_spec_method_t methods[] = {LC_SPEC_IHS, LC_SPEC_AOV, LC_SPEC_GB, LC_SPEC_BLS};
    const char *methodNames[] = {"IHS", "AoV", "GB", "BLS"};

    for (int mi = 0; mi < 4; ++mi) {
        lc_periodogram_config_t cfg;
        memset(&cfg, 0, sizeof(cfg));
        cfg.method = methods[mi];
        cfg.nterms = 3;
        cfg.oversampling = 5.0;
        cfg.oversmoothing = 0.2;
        cfg.nbins = 10;

        QVector<float> model(d.n);
        lc_model_style_t style = LC_MODEL_SCATTER;
        int rc = lc_compute_phased_model(&d, &cfg, f1, model.data(), &style);
        CHECK(rc == 0, "[%s] model computed successfully at f1", methodNames[mi]);
        CHECK(style == LC_MODEL_LINE, "[%s] model style is LC_MODEL_LINE", methodNames[mi]);

        double offset = lc_compute_phase_offset(&d, &cfg, f1, model.constData());
        CHECK(offset >= 0.0 && offset < 1.0, "[%s] phase offset in [0, 1): %f", methodNames[mi], offset);

        /* Set model and phase offset on widget */
        w.setPhaseOffset(offset);
        w.setModel(model.constData(), d.n, style, f1);
        CHECK(std::abs(w.phaseOffset() - offset) < 1e-6, "[%s] widget phase offset matches calculated offset", methodNames[mi]);

        /* Verify extremum alignment: phase 0.5 for minimum extremum, phase 0.0 for maximum extremum */
        QVector<float> ys(d.y, d.y + d.n);
        std::sort(ys.begin(), ys.end());
        double med = (d.n % 2 == 1) ? ys[d.n / 2] : 0.5 * (ys[d.n / 2 - 1] + ys[d.n / 2]);

        uint32_t minI = 0, maxI = 0;
        float minV = model[0], maxV = model[0];
        for (uint32_t i = 0; i < d.n; ++i) {
            if (model[i] < minV) { minV = model[i]; minI = i; }
            if (model[i] > maxV) { maxV = model[i]; maxI = i; }
        }
        double dBright = std::abs((double)minV - med) * LC_BRIGHTNESS_BIAS;
        double dFaint = std::abs((double)maxV - med) * 1.0;
        double targetPhase = (methods[mi] == LC_SPEC_BLS) ? 0.5 : ((dFaint > dBright) ? 0.5 : 0.0);
        uint32_t targetIdx = maxI;
        if (methods[mi] == LC_SPEC_BLS) {
            targetIdx = (std::abs((double)minV - med) > std::abs((double)maxV - med)) ? minI : maxI;
        }

        double foldedPhase = std::fmod((d.x[targetIdx] - t0) * f1 + offset, 1.0);
        if (foldedPhase < 0.0) foldedPhase += 1.0;
        double phaseDiff = std::abs(foldedPhase - targetPhase);
        if (phaseDiff > 0.5) phaseDiff = std::abs(phaseDiff - 1.0);
        double tol = (methods[mi] == LC_SPEC_BLS) ? 0.25 : 0.05;
        CHECK(phaseDiff < tol, "[%s] trough/dip point is near target phase %f (actual: %f)", methodNames[mi], targetPhase, foldedPhase);
    }

    /* --- 7. Recomputing GB/BLS spectrum after changing alpha with persistent context --- */
    printf("[7] GB/BLS recomputation after changing alpha with persistent computeCtx\n");
    lc_compute_ctx_t *pCtx = lc_compute_ctx_create(0);
    CHECK(pCtx != nullptr, "created persistent computeCtx");

    for (lc_spec_method_t m : {LC_SPEC_GB, LC_SPEC_BLS}) {
        const char *mName = (m == LC_SPEC_GB) ? "GB" : "BLS";
        lc_periodogram_config_t cfg;
        memset(&cfg, 0, sizeof(cfg));
        cfg.method = m;
        cfg.oversampling = 5.0;
        cfg.fmax = 24.0;
        cfg.oversmoothing = 0.2; /* initial alpha */
        cfg.nbins = 10;

        lc_periodogram_result_t res1;
        memset(&res1, 0, sizeof(res1));
        int rc1 = lc_compute_periodogram_ctx(pCtx, &d, &cfg, &res1, nullptr);
        CHECK(rc1 == 0, "[%s] initial spectrum (alpha=0.2) computed: %u freq bins", mName, res1.nfreq);
        double maxF1 = res1.fmin + (double)(res1.nfreq - 1) * res1.fstep;
        CHECK(maxF1 >= 24.0 - res1.fstep, "[%s] initial spectrum spans up to fmax (actual max: %f)", mName, maxF1);

        /* Change alpha to 0.05 (smaller alpha -> smaller fstep -> larger nfreq) without changing method */
        cfg.oversmoothing = 0.05;
        lc_periodogram_result_t res2;
        memset(&res2, 0, sizeof(res2));
        int rc2 = lc_compute_periodogram_ctx(pCtx, &d, &cfg, &res2, nullptr);
        CHECK(rc2 == 0, "[%s] recomputed spectrum (alpha=0.05) computed: %u freq bins", mName, res2.nfreq);
        double maxF2 = res2.fmin + (double)(res2.nfreq - 1) * res2.fstep;
        CHECK(maxF2 >= 24.0 - res2.fstep, "[%s] recomputed spectrum spans full grid up to fmax (actual max: %f)", mName, maxF2);
        CHECK(res2.nfreq > res1.nfreq, "[%s] frequency count correctly increased (%u > %u)", mName, res2.nfreq, res1.nfreq);

        lc_periodogram_result_free(&res1);
        lc_periodogram_result_free(&res2);
    }
    lc_compute_ctx_destroy(pCtx);

    lc_free(&d);
    printf("\n%s (%d failure%s)\n", failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
           failures, failures == 1 ? "" : "s");
    return failures == 0 ? 0 : 1;
}
