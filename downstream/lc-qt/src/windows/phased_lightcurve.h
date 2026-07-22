#pragma once

#include <QFrame>
#include <QPainter>
#include <QVector>
#include <QString>
#include <cmath>

extern "C" {
#include "lc_readout.h"
}

/*
 * PhasedLightCurveWidget — scatter plot of the light curve folded at the
 * pivot frequency.
 *
 * Phase is defined as (x[i] - x[0]) * freq mod 1 (time relative to the first
 * measurement, not the HJD value). The x-axis spans [0, 2]; every point is
 * drawn twice (at its phase and at phase + 1) so the wrap-around is visible.
 * The y-axis (magnitude, inverted: brighter on top) and overall formatting
 * match LightCurvePlotWidget exactly.
 *
 * Call setFrequency() whenever the spectrum pivot moves to re-fold the curve.
 */
class PhasedLightCurveWidget : public QFrame {
    Q_OBJECT
public:
    explicit PhasedLightCurveWidget(const QString &title, QWidget *parent = nullptr);
    ~PhasedLightCurveWidget() override;

    void setData(const lc_data_t *data);
    void clearData();
    bool hasData() const { return m_n > 0; }

public slots:
    /* Fold the curve at freq (cycles/day); freq <= 0 clears the folding. */
    void setFrequency(double freq);

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    QString m_title;
    QVector<double> m_x;
    QVector<float> m_y;
    unsigned int m_n = 0;
    double m_freq = -1.0;
};
