#pragma once

#include <QFrame>
#include <QPainter>
#include <QPoint>
#include <QElapsedTimer>
#include <QVector>
#include <QString>
#include <cmath>

/*
 * ZoomedSpectrumWidget — zoomed view of the NLL periodogram spectrum.
 *
 * Displays a sub-range (the FOV) of the full spectrum. A separate "pivot"
 * frequency (the indicator) acts as the zoom anchor:
 *  - the FOV is always kept inside the computed frequency range [fmin, fmax];
 *  - zooming (wheel or +/- buttons) re-centres the FOV on the pivot; when the
 *    pivot sits near an edge the FOV pins to that edge and the pivot appears
 *    off-centre rather than the view running past the computed range;
 *  - the pivot can be moved off-centre by binding it to the mouse (pressing
 *    LMB+RMB together toggles the binding) and moving the mouse.
 *
 * The display buffer is regenerated on every view change (zoom / pan / resize)
 * to maintain 1-8 points per pixel:
 *  - too dense (> 8 pt/px): peak-preserving min+max binning;
 *  - too sparse (< 0.75 pt/px): not-a-knot cubic spline at ~6 samples/px;
 *  - otherwise: direct copy.
 *
 * Interaction:
 *  - left click: move the pivot to the clicked frequency and centre the FOV on it;
 *  - left drag: translucent rubber-band rectangle; on release the pivot moves to
 *    the rectangle's centre and the FOV adopts its frequency range (the NLL axis
 *    stays auto-scaled);
 *  - LMB+RMB together: toggle binding the pivot to the mouse (while bound the
 *    pivot follows the pointer within the fixed FOV; the binding persists after
 *    the buttons are released);
 *  - mouse wheel: zoom by sqrt(2) per notch, re-centred on the pivot (the FOV
 *    never leaves the computed range);
 *  - zoomIn()/zoomOut() slots: zoom by factor 2 about the pivot;
 *  - stepFrequency()/multiplyFrequency() slots: move the pivot (arrow buttons,
 *    x2//2 and the "Other" multiplier menu). Both may move the pivot outside the
 *    computed range (the FOV itself never leaves it); stepFrequency() only keeps
 *    the frequency positive.
 */
class ZoomedSpectrumWidget : public QFrame {
    Q_OBJECT
public:
    explicit ZoomedSpectrumWidget(const QString &title, QWidget *parent = nullptr);
    ~ZoomedSpectrumWidget() override;

    void setFullSpectrum(double fmin, double fstep, const QVector<float> &nll);
    void clearData();
    bool hasData() const { return m_nfreq > 0; }

    /* Config zoom factor: applied when resetZoom == true (full-spectrum click). */
    void setZoomFactor(double factor) { m_zoomFactor = factor < 1.0 ? 1.0 : factor; }

    /* Move the pivot to freq and centre the FOV on it. resetZoom=true also resets
       zoom to m_zoomFactor. */
    void selectFrequency(double freq, bool resetZoom);

    /* Current view edges (frequency). */
    double viewMin() const;
    double viewMax() const;
    double centerFrequency() const { return m_pivotFreq; }

public slots:
    void zoomIn();  /* factor 2 */
    void zoomOut(); /* factor 2 */
    /* Shift the pivot by deltaFreq (arrow buttons); the FOV re-centres on it. */
    void stepFrequency(double deltaFreq);
    /* Scale the pivot frequency by factor (x2, /2, "Other" menu); FOV re-centres.
       Not clamped to the computed range: the pivot may leave it. */
    void multiplyFrequency(double factor);
    /* Adopt a rubber-band selection: the pivot moves to the selection centre and
       the FOV width becomes the selected frequency range, centred on the pivot.
       The NLL axis is left to auto-scale. */
    void setViewFromSelection(double fMin, double fMax);

signals:
    void centerFrequencyChanged(double freq);
    /* Middle mouse button click: peak-find at the clicked frequency within the current FOV. */
    void middleClicked(double freq, double fov);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private:
    QString m_title;

    /* Full spectrum data (copy). */
    double m_fmin = 0.0;
    double m_fstep = 1.0;
    QVector<float> m_nll;
    uint32_t m_nfreq = 0;

    /* View state. */
    double m_pivotFreq = -1.0; /* pivot (indicator) frequency; < 0 = no selection */
    double m_viewMin = 0.0;    /* FOV left edge (frequency); kept inside [fmin, fmax] */
    double m_zoomLevel = 1.0;  /* 1 = full range visible */
    double m_zoomFactor = 4.0; /* config: initial zoom on full-spectrum click */
    bool m_pivotBound = false; /* pivot bound to (tracking) the mouse */

    /* Rubber-band selection state. */
    QPoint m_pressPos;
    bool m_hasPress = false;
    bool m_dragging = false;
    QPoint m_dragCurrent;
    QElapsedTimer m_wheelIdleTimer; /* time since last wheel event (selection delay) */

    /* Regenerated display buffer (freq, value pairs at render density). */
    QVector<double> m_bufFreq;
    QVector<double> m_bufVal;

    double specFmin() const { return m_fmin; }
    double specFmax() const { return m_fmin + (double)(m_nfreq - 1) * m_fstep; }
    double freqAt(uint32_t i) const { return m_fmin + (double)i * m_fstep; }

    double fov() const; /* FOV width in frequency units */
    void centerOnPivot();
    void clampView();
    void regenerateBuffer();
    void viewChanged(); /* clamp + regen + emit + repaint */

    /* Pixel <-> data conversions for the current view (used by mouse handling). */
    double freqAtPixel(int px) const;
    void applyRubberBand(const QPoint &a, const QPoint &b);
    bool wheelIdle() const; /* true if no recent wheel event blocks the selection */
};
