#include "spectrum_plot.h"

#include <QPaintEvent>
#include <QPen>
#include <QFont>
#include <algorithm>
#include <cmath>

namespace {

constexpr double OVERZOOM_POINTS_PER_PIXEL = 0.75;
constexpr double SPLINE_SAMPLES_PER_PIXEL = 6.0;
constexpr int SPLINE_MARGIN_POINTS = 6;
constexpr int POINTS_PER_PIXEL_BUDGET = 8;

/* Not-a-knot cubic spline second derivatives M for strictly increasing xs (N >= 3). */
QVector<double> splineMoments(const QVector<double> &xs, const QVector<double> &ys) {
    const int N = xs.size();
    QVector<double> M(N, 0.0);
    if (N < 3) return M;
    const int n = N - 1; /* intervals */
    QVector<double> h(n);
    for (int i = 0; i < n; ++i) h[i] = xs[i + 1] - xs[i];

    const int m = n - 1; /* interior unknowns M[1..n-1] */
    if (m <= 0) return M;
    QVector<double> a(m), b(m), c(m), d(m);
    for (int i = 1; i <= n - 1; ++i) {
        int k = i - 1;
        double mu = h[i - 1] / (h[i - 1] + h[i]);
        double lam = 1.0 - mu;
        a[k] = mu;
        b[k] = 2.0;
        c[k] = lam;
        d[k] = 6.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]) / (h[i - 1] + h[i]);
    }
    /* Not-a-knot start: eliminate M[0] into row i=1. */
    {
        double h0 = h[0], h1 = h[1];
        b[0] += a[0] * (h0 + h1) / h1;
        c[0] -= a[0] * h0 / h1;
        a[0] = 0.0;
    }
    /* Not-a-knot end: eliminate M[n] into row i=n-1. */
    {
        double hn2 = h[n - 2], hn1 = h[n - 1];
        b[m - 1] += c[m - 1] * (hn2 + hn1) / hn2;
        a[m - 1] -= c[m - 1] * hn1 / hn2;
        c[m - 1] = 0.0;
    }
    /* Thomas algorithm. */
    for (int k = 1; k < m; ++k) {
        double w = a[k] / b[k - 1];
        b[k] -= w * c[k - 1];
        d[k] -= w * d[k - 1];
    }
    QVector<double> Mi(m);
    Mi[m - 1] = d[m - 1] / b[m - 1];
    for (int k = m - 2; k >= 0; --k) Mi[k] = (d[k] - c[k] * Mi[k + 1]) / b[k];
    for (int k = 0; k < m; ++k) M[k + 1] = Mi[k];
    M[0] = ((h[0] + h[1]) * M[1] - h[0] * M[2]) / h[1];
    M[n] = ((h[n - 2] + h[n - 1]) * M[n - 1] - h[n - 1] * M[n - 2]) / h[n - 2];
    return M;
}

double splineEval(const QVector<double> &xs, const QVector<double> &ys, const QVector<double> &M, double xp, int &seg) {
    const int N = xs.size();
    if (xp <= xs[0]) return ys[0];
    if (xp >= xs[N - 1]) return ys[N - 1];
    while (seg < N - 2 && xp > xs[seg + 1]) ++seg;
    while (seg > 0 && xp < xs[seg]) --seg;
    double h = xs[seg + 1] - xs[seg];
    if (h <= 0.0) return ys[seg];
    double A = (xs[seg + 1] - xp) / h;
    double B = (xp - xs[seg]) / h;
    return A * ys[seg] + B * ys[seg + 1] + ((A * A * A - A) * M[seg] + (B * B * B - B) * M[seg + 1]) * (h * h) / 6.0;
}

} // namespace

SpectrumPlotWidget::SpectrumPlotWidget(const QString &title, QWidget *parent) : QFrame(parent), m_title(title) {
    setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
    setMinimumHeight(150);
    setStyleSheet("background-color: #1a1a1f; border: 1px solid #32323d; border-radius: 4px;");
}

SpectrumPlotWidget::~SpectrumPlotWidget() = default;

void SpectrumPlotWidget::setProgressText(const QString &text) {
    if (m_progressText != text) {
        m_progressText = text;
        update();
    }
}

void SpectrumPlotWidget::setData(double fmin, double fstep, const QVector<float> &nll) {
    m_fmin = fmin;
    m_fstep = fstep > 0.0 ? fstep : 1.0;
    m_nll = nll;
    m_nfreq = (uint32_t)nll.size();
    update();
}

void SpectrumPlotWidget::clearData() {
    m_nll.clear();
    m_nfreq = 0;
    update();
}

void SpectrumPlotWidget::paintEvent(QPaintEvent *event) {
    QFrame::paintEvent(event);
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    int w = width();
    int h = height();
    int padL = 55, padR = 15, padT = 28, padB = 35;
    int plotW = w - padL - padR;
    int plotH = h - padT - padB;
    if (plotW <= 0 || plotH <= 0) return;

    // Title
    painter.setPen(QColor(203, 213, 225));
    painter.setFont(QFont("sans-serif", 9, QFont::Bold));
    painter.drawText(QRect(padL, 4, plotW, 20), Qt::AlignCenter, m_title);

    if (m_nfreq == 0) {
        painter.setPen(QPen(QColor(42, 42, 53), 1, Qt::DashLine));
        for (int i = 1; i < 6; i++) {
            int y = padT + i * plotH / 6;
            painter.drawLine(padL, y, padL + plotW, y);
            int x = padL + i * plotW / 6;
            painter.drawLine(x, padT, x, padT + plotH);
        }
        painter.setPen(QColor(100, 116, 139));
        painter.setFont(QFont("sans-serif", 11));
        painter.drawText(QRect(padL, padT, plotW, plotH), Qt::AlignCenter, "No data");
        // Progress overlay (show even when no spectrum yet)
        if (!m_progressText.isEmpty()) {
            painter.setPen(QColor(203, 213, 225));
            painter.setFont(QFont("sans-serif", 12, QFont::Bold));
            painter.drawText(QRect(padL, padT, plotW - 6, 22), Qt::AlignRight | Qt::AlignVCenter, m_progressText);
        }
        return;
    }

    // Data range — NLL y-axis starts at 0, 5% headroom above max
    double xMin = freqAt(0);
    double xMax = freqAt(m_nfreq - 1);
    float yMax = m_nll[0];
    for (uint32_t i = 1; i < m_nfreq; ++i) {
        if (m_nll[i] > yMax) yMax = m_nll[i];
    }
    if (!std::isfinite(yMax) || yMax <= 0.0f) yMax = 1.0f;
    float yMin = 0.0f;
    yMax += yMax * 0.05f; /* 5% margin above */
    float ySpan = yMax - yMin;

    double xSpan = xMax - xMin;
    if (xSpan < 1e-12) xSpan = 1.0;

    auto toPx = [&](double freq, float val) {
        double px = padL + ((freq - xMin) / xSpan) * plotW;
        double py = padT + plotH - ((double)(val - yMin) / ySpan) * plotH;
        return QPointF(px, py);
    };

    // Grid
    painter.setPen(QPen(QColor(42, 42, 53), 1, Qt::DashLine));
    const int gridN = 6;
    for (int i = 0; i <= gridN; i++) {
        int gy = padT + i * plotH / gridN;
        painter.drawLine(padL, gy, padL + plotW, gy);
        int gx = padL + i * plotW / gridN;
        painter.drawLine(gx, padT, gx, padT + plotH);
    }

    // Axes
    painter.setPen(QPen(QColor(100, 116, 139), 1.2));
    painter.drawLine(padL, padT + plotH, padL + plotW, padT + plotH);
    painter.drawLine(padL, padT, padL, padT + plotH);

    // Tick labels
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont("sans-serif", 9));
    int yDecimals = (int)std::floor(-std::log10((double)ySpan / 100.0));
    if (yDecimals < 0) yDecimals = 0;
    if (yDecimals > 6) yDecimals = 6;
    for (int i = 0; i <= gridN; i++) {
        float val = yMax - (float)i * ySpan / gridN;
        int gy = padT + i * plotH / gridN;
        painter.drawText(QRect(2, gy - 8, padL - 6, 16), Qt::AlignRight | Qt::AlignVCenter, QString::number(val, 'f', yDecimals));
    }
    int xDecimals = (std::fabs(xMax) > 1000.0 || std::fabs(xMin) > 1000.0) ? 0 : 3;
    for (int i = 0; i <= gridN; i += 2) {
        double val = xMin + (double)i * xSpan / gridN;
        int gx = padL + i * plotW / gridN;
        painter.drawText(QRect(gx - 30, padT + plotH + 4, 60, 16), Qt::AlignCenter, QString::number(val, 'f', xDecimals));
    }
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont("sans-serif", 8, QFont::Bold));
    painter.drawText(QRect(padL, h - 14, plotW, 14), Qt::AlignCenter, "Frequency [1/d]");

    // Build the rendered polyline
    QVector<QPointF> poly;
    double pointsPerPixel = (double)m_nfreq / (double)plotW;

    if (pointsPerPixel < OVERZOOM_POINTS_PER_PIXEL && m_nfreq >= 2) {
        // Sparse: local not-a-knot cubic spline over the (full) range + margin.
        int i0 = 0;
        int i1 = (int)m_nfreq;
        int j0 = std::max(0, i0 - SPLINE_MARGIN_POINTS);
        int j1 = std::min((int)m_nfreq, i1 + SPLINE_MARGIN_POINTS);
        int cnt = j1 - j0;
        if (cnt >= 3) {
            QVector<double> xs(cnt), ys(cnt);
            for (int k = 0; k < cnt; ++k) {
                xs[k] = freqAt((uint32_t)(j0 + k));
                ys[k] = (double)m_nll[j0 + k];
            }
            QVector<double> M = splineMoments(xs, ys);
            double a = std::max(xMin, xs[0]);
            double b = std::min(xMax, xs[cnt - 1]);
            int nSamples = std::max(64, (int)(plotW * SPLINE_SAMPLES_PER_PIXEL));
            poly.reserve(nSamples);
            int seg = 0;
            for (int s = 0; s < nSamples; ++s) {
                double xp = a + (b - a) * (double)s / (double)(nSamples - 1);
                double yp = splineEval(xs, ys, M, xp, seg);
                if (!std::isfinite(yp)) continue;
                if (yp < 0.0) yp = 0.0;
                poly.append(toPx(xp, (float)yp));
            }
        } else {
            for (uint32_t i = 0; i < m_nfreq; ++i) poly.append(toPx(freqAt(i), m_nll[i]));
        }
    } else if (pointsPerPixel > (double)POINTS_PER_PIXEL_BUDGET) {
        // Dense: peak-preserving downsample to <= 8 pts/pixel (min+max per bucket).
        int totalBuckets = (POINTS_PER_PIXEL_BUDGET / 2) * plotW; /* 2 pts (min+max) per bucket */
        if (totalBuckets < 1) totalBuckets = 1;
        double bucketSize = (double)m_nfreq / (double)totalBuckets;
        poly.reserve(totalBuckets * 2 + 2);
        for (int bkt = 0; bkt < totalBuckets; ++bkt) {
            int i0 = (int)((double)bkt * bucketSize);
            int i1 = (int)((double)(bkt + 1) * bucketSize);
            if (i1 > (int)m_nfreq) i1 = (int)m_nfreq;
            if (i0 >= i1) continue;
            int imin = i0, imax = i0;
            for (int i = i0; i < i1; ++i) {
                if (m_nll[i] < m_nll[imin]) imin = i;
                if (m_nll[i] > m_nll[imax]) imax = i;
            }
            if (imin <= imax) {
                poly.append(toPx(freqAt((uint32_t)imin), m_nll[imin]));
                if (imax != imin) poly.append(toPx(freqAt((uint32_t)imax), m_nll[imax]));
            } else {
                poly.append(toPx(freqAt((uint32_t)imax), m_nll[imax]));
                poly.append(toPx(freqAt((uint32_t)imin), m_nll[imin]));
            }
        }
    } else {
        // Moderate density: draw every point.
        poly.reserve((int)m_nfreq);
        for (uint32_t i = 0; i < m_nfreq; ++i) poly.append(toPx(freqAt(i), m_nll[i]));
    }

    // Draw the spectrum as a line
    if (poly.size() >= 2) {
        painter.setPen(QPen(QColor(59, 130, 246), 1.2)); /* #3b82f6 accent */
        painter.setBrush(Qt::NoBrush);
        painter.drawPolyline(poly.constData(), poly.size());
    } else if (poly.size() == 1) {
        painter.setPen(Qt::NoPen);
        painter.setBrush(QColor(59, 130, 246));
        painter.drawEllipse(poly[0], 1.5, 1.5);
    }

    // Progress overlay — upper-right corner of the plot area
    if (!m_progressText.isEmpty()) {
        painter.setPen(QColor(203, 213, 225));
        painter.setFont(QFont("sans-serif", 12, QFont::Bold));
        painter.drawText(QRect(padL, padT, plotW - 6, 22), Qt::AlignRight | Qt::AlignVCenter, m_progressText);
    }
}
