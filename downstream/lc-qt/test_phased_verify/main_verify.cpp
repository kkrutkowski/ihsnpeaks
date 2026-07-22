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
    const double t0 = d.x[0];
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

    lc_free(&d);
    printf("\n%s (%d failure%s)\n", failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
           failures, failures == 1 ? "" : "s");
    return failures == 0 ? 0 : 1;
}
