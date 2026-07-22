#include "zoomed_spectrum.h"
#include "spectrum_plot.h" /* SpectrumPlotWidget::PAD_* constants */

#include <QPaintEvent>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QResizeEvent>
#include <QPen>
#include <QFont>
#include <algorithm>
#include <cmath>

namespace {

constexpr double OVERZOOM_POINTS_PER_PIXEL = 0.75;
constexpr double SPLINE_SAMPLES_PER_PIXEL = 6.0;
constexpr int SPLINE_MARGIN_POINTS = 6;
constexpr int POINTS_PER_PIXEL_BUDGET = 8;
constexpr int DRAG_THRESHOLD_PX = 5;      /* click vs. rubber-band drag discrimination */
constexpr int RUBBERBAND_DELAY_MS = 500;  /* selection rectangle appears only after the wheel is idle this long */

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

ZoomedSpectrumWidget::ZoomedSpectrumWidget(const QString &title, QWidget *parent) : QFrame(parent), m_title(title) {
    setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
    setMinimumHeight(150);
    setStyleSheet("background-color: #1a1a1f; border: 1px solid #32323d; border-radius: 4px;");
}

ZoomedSpectrumWidget::~ZoomedSpectrumWidget() = default;

void ZoomedSpectrumWidget::setFullSpectrum(double fmin, double fstep, const QVector<float> &nll) {
    m_fmin = fmin;
    m_fstep = fstep > 0.0 ? fstep : 1.0;
    m_nll = nll;
    m_nfreq = (uint32_t)nll.size();
    if (m_nfreq > 0 && m_pivotFreq >= 0.0) {
        /* Keep existing pivot across recomputes: re-clamp and regenerate. */
        viewChanged();
    } else {
        m_bufFreq.clear();
        m_bufVal.clear();
        update();
    }
}

void ZoomedSpectrumWidget::clearData() {
    m_nll.clear();
    m_nfreq = 0;
    m_pivotFreq = -1.0;
    m_viewMin = 0.0;
    m_zoomLevel = 1.0;
    m_pivotBound = false;
    setMouseTracking(false);
    m_bufFreq.clear();
    m_bufVal.clear();
    update();
}

void ZoomedSpectrumWidget::selectFrequency(double freq, bool resetZoom) {
    if (m_nfreq == 0) return;
    m_pivotFreq = freq;
    if (resetZoom) m_zoomLevel = m_zoomFactor;
    centerOnPivot();
    viewChanged();
}

void ZoomedSpectrumWidget::zoomIn() {
    if (m_nfreq == 0 || m_pivotFreq < 0.0) return;
    m_zoomLevel *= 2.0;
    centerOnPivot();
    viewChanged();
}

void ZoomedSpectrumWidget::zoomOut() {
    if (m_nfreq == 0 || m_pivotFreq < 0.0) return;
    m_zoomLevel /= 2.0;
    centerOnPivot();
    viewChanged();
}

void ZoomedSpectrumWidget::stepFrequency(double deltaFreq) {
    /* Arrow buttons: shift the pivot (clamped to the computed range); the FOV
       re-centres on it via selectFrequency(). */
    if (m_nfreq == 0 || m_pivotFreq < 0.0) return;
    selectFrequency(std::clamp(m_pivotFreq + deltaFreq, specFmin(), specFmax()), false);
}

void ZoomedSpectrumWidget::multiplyFrequency(double factor) {
    /* x2, /2 and the "Other" menu: scale the pivot frequency. Unlike zooming,
       the multiplier buttons (and manual period edits) are deliberately NOT
       clamped to the computed range: the pivot may leave it, while clampView()
       still keeps the FOV itself inside [fmin, fmax]. */
    if (m_nfreq == 0 || m_pivotFreq < 0.0 || !(factor > 0.0)) return;
    selectFrequency(m_pivotFreq * factor, false);
}

double ZoomedSpectrumWidget::fov() const {
    double fullSpan = specFmax() - specFmin();
    return fullSpan / m_zoomLevel;
}

void ZoomedSpectrumWidget::centerOnPivot() {
    m_viewMin = m_pivotFreq - fov() / 2.0;
}

void ZoomedSpectrumWidget::clampView() {
    if (m_nfreq < 2) return;
    double lo = specFmin(), hi = specFmax();
    double fullSpan = hi - lo;
    if (fullSpan < 1e-12) {
        m_zoomLevel = 1.0;
        return;
    }

    /* Clamp zoom level to [1, fullSpan / minSpan]. */
    double minSpan = 8.0 * m_fstep;
    double maxZoom = fullSpan / minSpan;
    if (maxZoom < 1.0) maxZoom = 1.0;
    if (m_zoomLevel < 1.0) m_zoomLevel = 1.0;
    if (m_zoomLevel > maxZoom) m_zoomLevel = maxZoom;

    /* Keep the FOV inside the computed range [lo, hi]. When the pivot is near an
       edge the FOV pins to that edge (leaving the pivot off-centre) instead of
       running past the range. */
    double f = fov();
    if (f >= fullSpan) {
        m_viewMin = lo;
    } else {
        if (m_viewMin < lo) m_viewMin = lo;
        if (m_viewMin + f > hi) m_viewMin = hi - f;
    }
}

double ZoomedSpectrumWidget::viewMin() const {
    return m_viewMin;
}

double ZoomedSpectrumWidget::viewMax() const {
    return m_viewMin + fov();
}

double ZoomedSpectrumWidget::freqAtPixel(int px) const {
    int plotW = width() - SpectrumPlotWidget::PAD_L - SpectrumPlotWidget::PAD_R;
    if (plotW <= 0) return m_pivotFreq;
    double frac = (double)(px - SpectrumPlotWidget::PAD_L) / (double)plotW;
    frac = std::clamp(frac, 0.0, 1.0);
    return viewMin() + frac * (viewMax() - viewMin());
}

bool ZoomedSpectrumWidget::wheelIdle() const {
    /* Selection is allowed only once the wheel has been idle long enough. */
    return !m_wheelIdleTimer.isValid() || m_wheelIdleTimer.elapsed() >= RUBBERBAND_DELAY_MS;
}

void ZoomedSpectrumWidget::applyRubberBand(const QPoint &a, const QPoint &b) {
    if (m_nfreq == 0) return;
    int x0 = std::min(a.x(), b.x()), x1 = std::max(a.x(), b.x());
    /* Only the frequency range is used; the NLL axis auto-scales. */
    double f0 = freqAtPixel(x0), f1 = freqAtPixel(x1);
    setViewFromSelection(f0, f1);
}

void ZoomedSpectrumWidget::setViewFromSelection(double fMin, double fMax) {
    if (m_nfreq == 0) return;
    double lo = specFmin(), hi = specFmax();
    double fullSpan = hi - lo;
    if (fullSpan < 1e-12) return;
    if (fMin > fMax) std::swap(fMin, fMax);
    fMin = std::clamp(fMin, lo, hi);
    fMax = std::clamp(fMax, lo, hi);
    double sel = fMax - fMin;
    if (sel < 1e-12) {
        /* Degenerate (zero-width) selection: treat as a pan to that frequency. */
        selectFrequency(fMin, false);
        return;
    }
    m_zoomLevel = fullSpan / sel; /* clampView() constrains to [1, maxZoom] */
    /* Move the pivot to the middle of the selection, then centre the FOV on it:
       the FOV ends up exactly covering the selected frequency range. */
    m_pivotFreq = (fMin + fMax) / 2.0;
    centerOnPivot();
    viewChanged();
}

void ZoomedSpectrumWidget::regenerateBuffer() {
    m_bufFreq.clear();
    m_bufVal.clear();
    if (m_nfreq == 0 || m_pivotFreq < 0.0) return;

    int plotW = width() - SpectrumPlotWidget::PAD_L - SpectrumPlotWidget::PAD_R;
    if (plotW <= 0) return;

    double lo = specFmin(), hi = specFmax();
    double fullSpan = hi - lo;
    if (fullSpan < 1e-12) return;

    /* FOV (already clamped to [lo, hi] by clampView()). */
    double vMin = viewMin();
    double vMax = viewMax();
    double viewXSpan = vMax - vMin;
    if (viewXSpan < 1e-12) return;

    /* Data portion inside the view (clipping kept for safety). */
    double dataMin = std::max(lo, vMin);
    double dataMax = std::min(hi, vMax);
    if (dataMin >= dataMax) return;

    int i0 = (int)std::ceil((dataMin - m_fmin) / m_fstep);
    int i1 = (int)std::floor((dataMax - m_fmin) / m_fstep);
    if (i0 < 0) i0 = 0;
    if (i1 > (int)m_nfreq - 1) i1 = (int)m_nfreq - 1;
    int count = i1 - i0 + 1;
    if (count <= 0) return;

    /* Density relative to the pixels the data actually occupies. */
    double dataPixelWidth = (double)plotW * ((dataMax - dataMin) / viewXSpan);
    if (dataPixelWidth < 1.0) dataPixelWidth = 1.0;
    double pointsPerPixel = (double)count / dataPixelWidth;

    if (pointsPerPixel > (double)POINTS_PER_PIXEL_BUDGET) {
        /* Dense: peak-preserving min+max binning to <= 8 pts/pixel. */
        int totalBuckets = (int)((double)(POINTS_PER_PIXEL_BUDGET / 2) * dataPixelWidth);
        if (totalBuckets < 1) totalBuckets = 1;
        double bucketSize = (double)count / (double)totalBuckets;
        m_bufFreq.reserve(totalBuckets * 2);
        m_bufVal.reserve(totalBuckets * 2);
        for (int bkt = 0; bkt < totalBuckets; ++bkt) {
            int b0 = i0 + (int)((double)bkt * bucketSize);
            int b1 = i0 + (int)((double)(bkt + 1) * bucketSize);
            if (b1 > i1 + 1) b1 = i1 + 1;
            if (b0 >= b1) continue;
            int imin = b0, imax = b0;
            for (int i = b0; i < b1; ++i) {
                if (m_nll[i] < m_nll[imin]) imin = i;
                if (m_nll[i] > m_nll[imax]) imax = i;
            }
            if (imin <= imax) {
                m_bufFreq.append(freqAt((uint32_t)imin));
                m_bufVal.append((double)m_nll[imin]);
                if (imax != imin) {
                    m_bufFreq.append(freqAt((uint32_t)imax));
                    m_bufVal.append((double)m_nll[imax]);
                }
            } else {
                m_bufFreq.append(freqAt((uint32_t)imax));
                m_bufVal.append((double)m_nll[imax]);
                m_bufFreq.append(freqAt((uint32_t)imin));
                m_bufVal.append((double)m_nll[imin]);
            }
        }
    } else if (pointsPerPixel < OVERZOOM_POINTS_PER_PIXEL && count >= 2) {
        /* Sparse: not-a-knot cubic spline at ~6 samples/px. */
        int j0 = std::max(0, i0 - SPLINE_MARGIN_POINTS);
        int j1 = std::min((int)m_nfreq, i1 + 1 + SPLINE_MARGIN_POINTS);
        int cnt = j1 - j0;
        if (cnt >= 3) {
            QVector<double> xs(cnt), ys(cnt);
            for (int k = 0; k < cnt; ++k) {
                xs[k] = freqAt((uint32_t)(j0 + k));
                ys[k] = (double)m_nll[j0 + k];
            }
            QVector<double> M = splineMoments(xs, ys);
            double a = std::max(dataMin, xs[0]);
            double b = std::min(dataMax, xs[cnt - 1]);
            int nSamples = std::max(64, (int)(dataPixelWidth * SPLINE_SAMPLES_PER_PIXEL));
            m_bufFreq.reserve(nSamples);
            m_bufVal.reserve(nSamples);
            int seg = 0;
            for (int s = 0; s < nSamples; ++s) {
                double xp = a + (b - a) * (double)s / (double)(nSamples - 1);
                double yp = splineEval(xs, ys, M, xp, seg);
                if (!std::isfinite(yp)) continue;
                if (yp < 0.0) yp = 0.0;
                m_bufFreq.append(xp);
                m_bufVal.append(yp);
            }
        } else {
            for (int i = i0; i <= i1; ++i) {
                m_bufFreq.append(freqAt((uint32_t)i));
                m_bufVal.append((double)m_nll[i]);
            }
        }
    } else {
        /* Moderate density: direct copy. */
        m_bufFreq.reserve(count);
        m_bufVal.reserve(count);
        for (int i = i0; i <= i1; ++i) {
            m_bufFreq.append(freqAt((uint32_t)i));
            m_bufVal.append((double)m_nll[i]);
        }
    }
}

void ZoomedSpectrumWidget::viewChanged() {
    clampView();
    regenerateBuffer();
    if (m_pivotFreq >= 0.0)
        emit centerFrequencyChanged(m_pivotFreq);
    update();
}

void ZoomedSpectrumWidget::mousePressEvent(QMouseEvent *event) {
    /* Pressing LMB+RMB together toggles binding the pivot to the mouse. The
       binding persists after the buttons are released: while bound, mouse
       tracking is enabled so move events arrive with no button held. */
    if ((event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::RightButton)) {
        if (m_nfreq > 0 && m_pivotFreq >= 0.0) {
            m_pivotBound = !m_pivotBound;
            setMouseTracking(m_pivotBound);
            m_hasPress = false;
            m_dragging = false;
            update();
        }
        return;
    }
    if (event->button() == Qt::LeftButton && m_nfreq > 0 && m_pivotFreq >= 0.0) {
        m_pressPos = event->position().toPoint();
        m_dragCurrent = m_pressPos;
        m_hasPress = true;
        m_dragging = false;
        return;
    }
    QFrame::mousePressEvent(event);
}

void ZoomedSpectrumWidget::mouseMoveEvent(QMouseEvent *event) {
    if (m_pivotBound && m_nfreq > 0 && m_pivotFreq >= 0.0) {
        /* Pivot tracks the pointer within the (fixed) FOV, so it can be moved
           off-centre. The FOV itself does not move while bound. */
        double lo = specFmin(), hi = specFmax();
        double f = std::clamp(freqAtPixel(event->position().x()), lo, hi);
        if (f != m_pivotFreq) {
            m_pivotFreq = f;
            emit centerFrequencyChanged(m_pivotFreq);
            update();
        }
        return;
    }
    if (m_hasPress) {
        /* The selection rectangle only appears once the wheel has been idle for
           RUBBERBAND_DELAY_MS. */
        if (!m_dragging && wheelIdle() && (event->position().toPoint() - m_pressPos).manhattanLength() > DRAG_THRESHOLD_PX)
            m_dragging = true;
        if (m_dragging) {
            m_dragCurrent = event->position().toPoint();
            update();
            return;
        }
    }
    QFrame::mouseMoveEvent(event);
}

void ZoomedSpectrumWidget::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton && m_hasPress) {
        bool wasDrag = m_dragging;
        QPoint releasePos = event->position().toPoint();
        m_hasPress = false;
        m_dragging = false;
        if (wasDrag) {
            applyRubberBand(m_pressPos, releasePos);
        } else {
            /* Plain click: move the pivot to the clicked frequency and centre the
               FOV on it (clamped to the data range). */
            double lo = specFmin(), hi = specFmax();
            selectFrequency(std::clamp(freqAtPixel(releasePos.x()), lo, hi), false);
        }
        update();
        return;
    }
    QFrame::mouseReleaseEvent(event);
}

void ZoomedSpectrumWidget::wheelEvent(QWheelEvent *event) {
    if (m_nfreq == 0 || m_pivotFreq < 0.0) {
        QFrame::wheelEvent(event);
        return;
    }
    int delta = event->angleDelta().y();
    if (delta == 0) {
        QFrame::wheelEvent(event);
        return;
    }

    /* Each scroll resets the marked (selection) area; the selection rectangle
       only reappears after the wheel has been idle for RUBBERBAND_DELAY_MS. */
    m_wheelIdleTimer.restart();
    if (m_dragging) m_dragging = false;

    double lo = specFmin(), hi = specFmax();
    double fullSpan = hi - lo;
    if (fullSpan < 1e-12) {
        QFrame::wheelEvent(event);
        return;
    }

    /* Pivot-anchored zoom: resize the FOV and re-centre it on the pivot. The zoom
       does NOT follow the cursor; the pivot only moves via click or mouse-binding.
       clampView() keeps the FOV inside [lo, hi] (pinning it to an edge when the
       pivot is near one, which leaves the pivot off-centre). */
    double f = fov();
    double newFov = f / std::pow(std::sqrt(2.0), (double)delta / 120.0);
    double minSpan = 8.0 * m_fstep;
    if (newFov < minSpan) newFov = minSpan;
    if (newFov > fullSpan) newFov = fullSpan;
    m_zoomLevel = fullSpan / newFov;
    centerOnPivot();
    viewChanged();
    event->accept();
}

void ZoomedSpectrumWidget::resizeEvent(QResizeEvent *event) {
    QFrame::resizeEvent(event);
    if (m_nfreq > 0 && m_pivotFreq >= 0.0) {
        regenerateBuffer();
        update();
    }
}

void ZoomedSpectrumWidget::paintEvent(QPaintEvent *event) {
    QFrame::paintEvent(event);
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    int w = width();
    int h = height();
    int padL = SpectrumPlotWidget::PAD_L, padR = SpectrumPlotWidget::PAD_R;
    int padT = SpectrumPlotWidget::PAD_T, padB = SpectrumPlotWidget::PAD_B;
    int plotW = w - padL - padR;
    int plotH = h - padT - padB;
    if (plotW <= 0 || plotH <= 0) return;

    // Title
    painter.setPen(QColor(203, 213, 225));
    painter.setFont(QFont("sans-serif", 9, QFont::Bold));
    painter.drawText(QRect(padL, 4, plotW, 20), Qt::AlignCenter, m_title);

    bool hasSelection = (m_nfreq > 0 && m_pivotFreq >= 0.0);

    if (!hasSelection) {
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
        return;
    }

    // View window (x range) — clamped to the computed range by clampView().
    double xMin = viewMin();
    double xMax = viewMax();
    double xSpan = xMax - xMin;
    if (xSpan < 1e-12) xSpan = 1.0;

    // Y range — always auto-scaled from 0 to the visible maximum + 5% headroom.
    double yMax = 1.0;
    if (m_bufVal.size() > 0) {
        yMax = m_bufVal[0];
        for (int i = 1; i < m_bufVal.size(); ++i)
            if (m_bufVal[i] > yMax) yMax = m_bufVal[i];
    }
    if (!std::isfinite(yMax) || yMax <= 0.0) yMax = 1.0;
    yMax += yMax * 0.05; /* 5% margin above */
    double yMin = 0.0;
    double ySpan = yMax - yMin;
    if (ySpan < 1e-12) ySpan = 1.0;

    auto toPx = [&](double freq, double val) {
        double px = padL + ((freq - xMin) / xSpan) * plotW;
        double py = padT + plotH - ((val - yMin) / ySpan) * plotH;
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
    int yDecimals = (int)std::floor(-std::log10(ySpan / 100.0));
    if (yDecimals < 0) yDecimals = 0;
    if (yDecimals > 6) yDecimals = 6;
    for (int i = 0; i <= gridN; i++) {
        double val = yMax - (double)i * ySpan / gridN;
        int gy = padT + i * plotH / gridN;
        painter.drawText(QRect(2, gy - 8, padL - 6, 16), Qt::AlignRight | Qt::AlignVCenter, QString::number(val, 'f', yDecimals));
    }
    int xDecimals = 3;
    if (std::fabs(xMax) > 1000.0 || std::fabs(xMin) > 1000.0) {
        xDecimals = 0;
    } else if (xSpan < 1.0) {
        /* Deep zoom: increase precision so tick labels remain distinguishable. */
        xDecimals = (int)std::ceil(-std::log10(xSpan)) + 2;
        if (xDecimals > 8) xDecimals = 8;
    }
    for (int i = 0; i <= gridN; i += 2) {
        double val = xMin + (double)i * xSpan / gridN;
        int gx = padL + i * plotW / gridN;
        painter.drawText(QRect(gx - 30, padT + plotH + 4, 60, 16), Qt::AlignCenter, QString::number(val, 'f', xDecimals));
    }
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont("sans-serif", 8, QFont::Bold));
    painter.drawText(QRect(padL, h - 14, plotW, 14), Qt::AlignCenter, "Frequency [1/d]");

    // Draw the spectrum from the pre-processed display buffer
    if (m_bufFreq.size() >= 2) {
        QVector<QPointF> poly(m_bufFreq.size());
        for (int i = 0; i < m_bufFreq.size(); ++i)
            poly[i] = toPx(m_bufFreq[i], m_bufVal[i]);
        painter.setPen(QPen(QColor(59, 130, 246), 1.2)); /* #3b82f6 accent */
        painter.setBrush(Qt::NoBrush);
        painter.drawPolyline(poly.constData(), poly.size());
    } else if (m_bufFreq.size() == 1) {
        painter.setPen(Qt::NoPen);
        painter.setBrush(QColor(59, 130, 246));
        painter.drawEllipse(toPx(m_bufFreq[0], m_bufVal[0]), 1.5, 1.5);
    }

    // Pivot indicator — thin vertical line (may be off-centre when the FOV is
    // pinned to an edge or the pivot has been moved with the mouse).
    if (m_pivotFreq >= xMin && m_pivotFreq <= xMax) {
        int ix = padL + (int)std::round(((m_pivotFreq - xMin) / xSpan) * plotW);
        painter.setPen(QPen(QColor(245, 158, 11), 1)); /* #f59e0b amber */
        painter.drawLine(ix, padT, ix, padT + plotH);
    }

    // Rubber-band selection overlay — translucent rectangle spanning the full
    // plot height (only the frequency range is meaningful).
    if (m_dragging) {
        int x0 = std::min(m_pressPos.x(), m_dragCurrent.x());
        int x1 = std::max(m_pressPos.x(), m_dragCurrent.x());
        QRect r(x0, padT, x1 - x0, plotH);
        painter.setPen(QPen(QColor(59, 130, 246, 200), 1));
        painter.setBrush(QColor(59, 130, 246, 50));
        painter.drawRect(r);
    }
}
