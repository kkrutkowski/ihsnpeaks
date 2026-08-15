#include "phased_lightcurve.h"

#include <QPaintEvent>
#include <QPen>
#include <QFont>
#include <QPainterPath>
#include <algorithm>
#include <cmath>
#include <vector>

PhasedLightCurveWidget::PhasedLightCurveWidget(const QString &title, QWidget *parent)
    : QFrame(parent), m_title(title) {
    setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
    setMinimumHeight(150);
    setStyleSheet("background-color: #1a1a1f; border: 1px solid #32323d; border-radius: 4px;");
}

PhasedLightCurveWidget::~PhasedLightCurveWidget() = default;

void PhasedLightCurveWidget::setData(const lc_data_t *data) {
    if (!data || data->n == 0) {
        clearData();
        return;
    }
    m_n = data->n;
    m_x.resize(m_n);
    m_y.resize(m_n);
    for (unsigned int i = 0; i < m_n; i++) {
        m_x[i] = data->x[i];
        m_y[i] = data->y[i];
    }
    m_phaseOffset = 0.0;
    update();
}

void PhasedLightCurveWidget::clearData() {
    m_x.clear();
    m_y.clear();
    m_n = 0;
    m_phaseOffset = 0.0;
    update();
}

void PhasedLightCurveWidget::setFrequency(double freq) {
    m_freq = freq;
    update();
}

void PhasedLightCurveWidget::setPhaseOffset(double offset) {
    double wrapped = std::fmod(offset, 1.0);
    if (wrapped < 0.0) wrapped += 1.0;
    if (m_phaseOffset != wrapped) {
        m_phaseOffset = wrapped;
        update();
    }
}

void PhasedLightCurveWidget::setModel(const float *model, unsigned int n, lc_model_style_t style, double freq) {
    if (!model || n == 0 || n != m_n) {
        clearModel();
        return;
    }
    m_model.resize(n);
    for (unsigned int i = 0; i < n; i++) m_model[i] = model[i];
    m_modelStyle = style;
    m_modelFreq = freq;
    m_hasModel = true;
    /* Coefficient of determination: fraction of the data variance explained by
       the model (R^2 = 1 - SS_res / SS_tot). */
    double mean = 0.0;
    for (unsigned int i = 0; i < n; i++) mean += (double)m_y[i];
    mean /= (double)n;
    double ssTot = 0.0, ssRes = 0.0;
    for (unsigned int i = 0; i < n; i++) {
        double dTot = (double)m_y[i] - mean;
        ssTot += dTot * dTot;
        double dRes = (double)m_y[i] - (double)m_model[i];
        ssRes += dRes * dRes;
    }
    m_r2 = (ssTot > 0.0) ? (1.0 - ssRes / ssTot) : 0.0;
    update();
}

void PhasedLightCurveWidget::clearModel() {
    m_model.clear();
    m_hasModel = false;
    m_modelFreq = -1.0;
    m_r2 = 0.0;
    update();
}

void PhasedLightCurveWidget::setDisplayModel(bool show) {
    if (m_displayModel != show) {
        m_displayModel = show;
        update();
    }
}

void PhasedLightCurveWidget::paintEvent(QPaintEvent *event) {
    QFrame::paintEvent(event);
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    int w = width();
    int h = height();
    int padL = 55, padR = 15, padT = 28, padB = 35;
    int plotW = w - padL - padR;
    int plotH = h - padT - padB;

    if (plotW <= 0 || plotH <= 0) return;

    // Draw title
    painter.setPen(QColor(203, 213, 225));
    painter.setFont(QFont("sans-serif", 9, QFont::Bold));
    painter.drawText(QRect(padL, 4, plotW, 20), Qt::AlignCenter, m_title);

    if (m_n == 0 || !(m_freq > 0.0)) {
        // Empty state: grid + "No data"
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

    // Compute magnitude range (same as the raw light curve widget)
    float yMin = m_y[0], yMax = m_y[0];
    for (unsigned int i = 1; i < m_n; i++) {
        if (m_y[i] < yMin) yMin = m_y[i];
        if (m_y[i] > yMax) yMax = m_y[i];
    }

    // Add 5% margin to y range
    float ySpan = yMax - yMin;
    if (ySpan < 1e-6f) ySpan = 1.0f;
    yMin -= ySpan * 0.05f;
    yMax += ySpan * 0.05f;
    ySpan = yMax - yMin;

    // Fixed phase range [0, 2]
    const double xMin = 0.0, xMax = 2.0, xSpan = 2.0;

    // Draw grid
    painter.setPen(QPen(QColor(42, 42, 53), 1, Qt::DashLine));
    const int gridN = 6;
    for (int i = 0; i <= gridN; i++) {
        int gy = padT + i * plotH / gridN;
        painter.drawLine(padL, gy, padL + plotW, gy);
        int gx = padL + i * plotW / gridN;
        painter.drawLine(gx, padT, gx, padT + plotH);
    }

    // Draw axes
    painter.setPen(QPen(QColor(100, 116, 139), 1.2));
    painter.drawLine(padL, padT + plotH, padL + plotW, padT + plotH); // X
    painter.drawLine(padL, padT, padL, padT + plotH);                 // Y

    // Axis tick labels
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont("sans-serif", 9));

    // Y-axis labels (inverted: brighter/smaller mag on top)
    int yDecimals = (int)std::floor(-std::log10(ySpan / 100.0));
    if (yDecimals < 0) yDecimals = 0;
    for (int i = 0; i <= gridN; i++) {
        float val = yMin + (float)i * ySpan / gridN;
        int gy = padT + i * plotH / gridN;
        painter.drawText(QRect(2, gy - 8, padL - 6, 16), Qt::AlignRight | Qt::AlignVCenter,
                         QString::number(val, 'f', yDecimals));
    }

    // X-axis labels: integer multiples of 0.25, always 2 decimals
    for (int k = 0; k <= 8; k++) {
        double val = 0.25 * k;
        int gx = padL + (int)std::round((val / xSpan) * plotW);
        painter.drawText(QRect(gx - 30, padT + plotH + 4, 60, 16), Qt::AlignCenter,
                         QString::number(val, 'f', 2));
    }

    // X-axis title
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont("sans-serif", 8, QFont::Bold));
    painter.drawText(QRect(padL, h - 14, plotW, 14), Qt::AlignCenter, "Phase");

    // Draw data points as scatter: each point at its phase and at phase + 1
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(59, 130, 246)); // #3b82f6 accent

    const double t0 = (m_x[0] + m_x[m_n - 1]) / 2.0;
    for (unsigned int i = 0; i < m_n; i++) {
        double phase = std::fmod((m_x[i] - t0) * m_freq + m_phaseOffset, 1.0);
        if (phase < 0.0) phase += 1.0;
        double py = padT + ((double)(m_y[i] - yMin) / ySpan) * plotH;
        double px = padL + ((phase - xMin) / xSpan) * plotW;
        painter.drawEllipse(QPointF(px, py), 1.5, 1.5);
        px = padL + ((phase + 1.0 - xMin) / xSpan) * plotW;
        painter.drawEllipse(QPointF(px, py), 1.5, 1.5);
    }

    // Draw model overlay (semi-transparent, distinct colour)
    if (m_displayModel && m_hasModel && m_model.size() == (int)m_n) {
        QColor modelColor(245, 158, 11, 200); // amber #f59e0b, alpha ~200

        if (m_modelStyle == LC_MODEL_SCATTER) {
            // Scatter: draw ellipses at model positions
            painter.setPen(Qt::NoPen);
            painter.setBrush(modelColor);
            for (unsigned int i = 0; i < m_n; i++) {
                double phase = std::fmod((m_x[i] - t0) * m_modelFreq + m_phaseOffset, 1.0);
                if (phase < 0.0) phase += 1.0;
                double py = padT + ((double)(m_model[i] - yMin) / ySpan) * plotH;
                double px = padL + ((phase - xMin) / xSpan) * plotW;
                painter.drawEllipse(QPointF(px, py), 1.5, 1.5);
                px = padL + ((phase + 1.0 - xMin) / xSpan) * plotW;
                painter.drawEllipse(QPointF(px, py), 1.5, 1.5);
            }
        } else {
            // Line: sort by phase, draw connected polyline
            struct PhasePoint {
                double phase;
                float val;
            };
            std::vector<PhasePoint> pts;
            pts.reserve(m_n);
            for (unsigned int i = 0; i < m_n; i++) {
                double phase = std::fmod((m_x[i] - t0) * m_modelFreq + m_phaseOffset, 1.0);
                if (phase < 0.0) phase += 1.0;
                pts.push_back({phase, m_model[i]});
            }
            std::sort(pts.begin(), pts.end(), [](const PhasePoint &a, const PhasePoint &b) {
                return a.phase < b.phase;
            });

            painter.setPen(QPen(modelColor, 3.0));
            painter.setBrush(Qt::NoBrush);

            // Draw the model curve for [0,1] and [1,2] (wrap)
            for (int wrap = 0; wrap < 2; wrap++) {
                QPainterPath path;
                bool first = true;
                for (const auto &pt : pts) {
                    double ph = pt.phase + (double)wrap;
                    double px = padL + ((ph - xMin) / xSpan) * plotW;
                    double py = padT + ((double)(pt.val - yMin) / ySpan) * plotH;
                    if (first) {
                        path.moveTo(px, py);
                        first = false;
                    } else {
                        path.lineTo(px, py);
                    }
                }
                painter.drawPath(path);
            }
        }
    }

    /* R^2 overlay — top-right, same style as the computation percentage. */
    if (m_displayModel && m_hasModel) {
        painter.setPen(QColor(203, 213, 225));
        painter.setFont(QFont("sans-serif", 12, QFont::Bold));
        painter.drawText(QRect(padL, padT, plotW - 6, 22), Qt::AlignRight | Qt::AlignVCenter,
                         QString("R\u00B2 = %1").arg(m_r2, 0, 'f', 3));
    }
}
