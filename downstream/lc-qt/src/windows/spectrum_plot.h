#pragma once

#include <QFrame>
#include <QPainter>
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

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    QString m_title;
    QString m_progressText;
    double m_fmin = 0.0;
    double m_fstep = 1.0;
    QVector<float> m_nll;
    uint32_t m_nfreq = 0;

    double freqAt(uint32_t i) const { return m_fmin + (double)i * m_fstep; }
};
