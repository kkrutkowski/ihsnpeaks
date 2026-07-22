#pragma once

#include <QFrame>
#include <QPainter>
#include <QPoint>
#include <QVector>
#include <QString>
#include <cmath>

/*
 * SpectrumPlotWidget — line plot of a periodogram (NLL vs frequency).
 *
 * Rendering follows src/debug/spec_viewer.py:
 *  - dense (>= ~0.75 pt/px): peak-preserving downsample to <= 8 pts/pixel
 *    (min+max per bucket so local maxima of neighboring frequencies survive);
 *  - sparse (< 0.75 pt/px, e.g. future zoom): not-a-knot cubic spline at ~6 samples/px.
 * The y (NLL) axis uses a 5% margin.
 *
 * Interaction:
 *  - left click: selects a frequency (emits frequencyClicked);
 *  - left drag: translucent rubber-band rectangle; on release emits
 *    rangeSelected with the selected frequency range (the NLL axis of the
 *    zoomed view auto-scales; only frequency is constrained).
 */
class SpectrumPlotWidget : public QFrame {
    Q_OBJECT
public:
    explicit SpectrumPlotWidget(const QString &title, QWidget *parent = nullptr);
    ~SpectrumPlotWidget() override;

    void setData(double fmin, double fstep, const QVector<float> &nll);
    void clearData();
    bool hasData() const { return m_nfreq > 0; }

    void setProgressText(const QString &text);

    /* Selected-frequency indicator (thin vertical line); freq < 0 clears it. */
    void setSelectedFrequency(double freq);
    void clearSelectedFrequency() { setSelectedFrequency(-1.0); }

    /* Plot-area padding constants (shared with ZoomedSpectrumWidget). */
    static constexpr int PAD_L = 55;
    static constexpr int PAD_R = 15;
    static constexpr int PAD_T = 28;
    static constexpr int PAD_B = 35;

signals:
    void frequencyClicked(double freq);
    /* Rubber-band selection: frequency range only (NLL axis auto-scales). */
    void rangeSelected(double fMin, double fMax);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;

private:
    QString m_title;
    QString m_progressText;
    double m_fmin = 0.0;
    double m_fstep = 1.0;
    QVector<float> m_nll;
    uint32_t m_nfreq = 0;
    double m_selectedFreq = -1.0;

    /* Rubber-band selection state. */
    QPoint m_pressPos;
    bool m_hasPress = false;
    bool m_dragging = false;
    QPoint m_dragCurrent;

    double freqAt(uint32_t i) const { return m_fmin + (double)i * m_fstep; }
    double freqAtPixel(int px) const;
};
