#ifndef PERIOD_WORKER_H
#define PERIOD_WORKER_H

#include <QObject>
#include <QString>
#include <QVector>
#include <QPair>

extern "C" {
#include "lc_readout.h"
#include "lc_period.h"
}

/* Deep copy of an lc_data_t for use on the worker thread (frees only x/y/dy; _buf stays null). */
lc_data_t cloneLcData(const lc_data_t *src);

/* Format coordinate (frequency or period) following upstream precision conventions. */
QString formatCoordinateUpstream(double freq, double deltaT, double oversampling, int nterms, bool formatFrequency);

/* Worker that runs the (blocking, multithreaded) periodogram computation off the GUI thread. */
class PeriodogramTask : public QObject {
    Q_OBJECT
public:
    PeriodogramTask(lc_compute_ctx_t *ctx, const lc_data_t *data, const lc_periodogram_config_t &cfg, lc_progress_t *progress);
    ~PeriodogramTask() override;

public slots:
    void run();

signals:
    void finished(double fmin, double fstep, QVector<float> nll, QVector<QPair<double, float>> peaks);
    void cancelled();
    void failed(QString msg);

private:
    lc_compute_ctx_t *m_ctx;
    lc_data_t m_data;
    lc_periodogram_config_t m_cfg;
    lc_progress_t *m_progress;
};

/* Worker that computes the phased model overlay off the GUI thread so the UI
   stays responsive while the pivot is being scrolled. */
class PhasedModelWorker : public QObject {
    Q_OBJECT
public:
    /* Runs on the worker thread; takes ownership of the cloned data. */
    void compute(lc_data_t data, lc_periodogram_config_t cfg, double freq);

signals:
    void done(int rc, double freq, QVector<float> model, int style, qint64 elapsedNs);
};

#endif // PERIOD_WORKER_H
