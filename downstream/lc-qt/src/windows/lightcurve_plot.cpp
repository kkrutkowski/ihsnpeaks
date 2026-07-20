#include "lightcurve_plot.h"

#include <QPaintEvent>
#include <QPen>
#include <QFont>
#include <algorithm>
#include <cmath>

LightCurvePlotWidget::LightCurvePlotWidget(const QString &title, QWidget *parent)
    : QFrame(parent), m_title(title) {
    setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
    setMinimumHeight(150);
    setStyleSheet("background-color: #1a1a1f; border: 1px solid #32323d; border-radius: 4px;");
}

LightCurvePlotWidget::~LightCurvePlotWidget() = default;

void LightCurvePlotWidget::setData(const lc_data_t *data) {
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
    update();
}

void LightCurvePlotWidget::clearData() {
    m_x.clear();
    m_y.clear();
    m_n = 0;
    update();
}

void LightCurvePlotWidget::paintEvent(QPaintEvent *event) {
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

    if (m_n == 0) {
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

    // Compute data range
    double xMin = m_x[0], xMax = m_x[0];
    float yMin = m_y[0], yMax = m_y[0];
    for (unsigned int i = 1; i < m_n; i++) {
        if (m_x[i] < xMin) xMin = m_x[i];
        if (m_x[i] > xMax) xMax = m_x[i];
        if (m_y[i] < yMin) yMin = m_y[i];
        if (m_y[i] > yMax) yMax = m_y[i];
    }

    // Add 5% margin to y range
    float ySpan = yMax - yMin;
    if (ySpan < 1e-6f) ySpan = 1.0f;
    yMin -= ySpan * 0.05f;
    yMax += ySpan * 0.05f;
    ySpan = yMax - yMin;

    double xSpan = xMax - xMin;
    if (xSpan < 1e-9) xSpan = 1.0;

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

    // Y-axis labels: precision = 1/100 of visible amplitude
    int yDecimals = (int)std::floor(-std::log10(ySpan / 100.0));
    if (yDecimals < 0) yDecimals = 0;
    for (int i = 0; i <= gridN; i++) {
        float val = yMax - (float)i * ySpan / gridN;
        int gy = padT + i * plotH / gridN;
        painter.drawText(QRect(2, gy - 8, padL - 6, 16), Qt::AlignRight | Qt::AlignVCenter,
                         QString::number(val, 'f', yDecimals));
    }

    // X-axis labels: up to 3 decimals, but if values > 1000 show integers only
    int xDecimals;
    if (std::fabs(xMax) > 1000.0 || std::fabs(xMin) > 1000.0) {
        xDecimals = 0;
    } else {
        xDecimals = 3;
    }
    for (int i = 0; i <= gridN; i += 2) {
        double val = xMin + (double)i * xSpan / gridN;
        int gx = padL + i * plotW / gridN;
        painter.drawText(QRect(gx - 30, padT + plotH + 4, 60, 16), Qt::AlignCenter,
                         QString::number(val, 'f', xDecimals));
    }

    // X-axis title
    painter.setPen(QColor(148, 163, 184));
    painter.setFont(QFont("sans-serif", 8, QFont::Bold));
    painter.drawText(QRect(padL, h - 14, plotW, 14), Qt::AlignCenter, "Time [d]");

    // Draw data points as scatter
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(59, 130, 246)); // #3b82f6 accent

    for (unsigned int i = 0; i < m_n; i++) {
        double px = padL + ((m_x[i] - xMin) / xSpan) * plotW;
        double py = padT + ((double)(yMax - m_y[i]) / ySpan) * plotH;
        painter.drawEllipse(QPointF(px, py), 1.5, 1.5);
    }
}
