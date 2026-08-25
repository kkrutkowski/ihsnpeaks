#include "period_worker.h"

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <QElapsedTimer>

lc_data_t cloneLcData(const lc_data_t *src) {
    lc_data_t copy;
    memset(&copy, 0, sizeof(copy));
    if (!src || src->n == 0) return copy;
    copy.n = src->n;
    copy.magnitude = src->magnitude;
    copy._buf = nullptr;
    copy.x = (double *)malloc((size_t)src->n * sizeof(double));
    copy.y = (float *)malloc((size_t)src->n * sizeof(float));
    copy.dy = (float *)malloc((size_t)src->n * sizeof(float));
    if (copy.x) memcpy(copy.x, src->x, (size_t)src->n * sizeof(double));
    if (copy.y) memcpy(copy.y, src->y, (size_t)src->n * sizeof(float));
    if (copy.dy) memcpy(copy.dy, src->dy, (size_t)src->n * sizeof(float));
    return copy;
}

QString formatCoordinateUpstream(double freq, double deltaT, double oversampling, int nterms, bool formatFrequency) {
    if (!(freq > 0.0) || !std::isfinite(freq)) return QStringLiteral("inf");
    double arg = deltaT * oversampling * (double)nterms;
    int precision = (arg > 0.0) ? 1 + (int)std::log10(arg) : 1;
    int sig = precision + (int)std::ceil(std::log10(freq));
    if (sig < 1) sig = 1;
    if (sig > 17) sig = 17;
    double val = formatFrequency ? freq : (1.0 / freq);
    int dec;
    if (val >= 1.0) {
        int intDigits = (int)std::ceil(std::log10(val));
        if (intDigits < 0) intDigits = 0;
        dec = sig - intDigits;
        if (dec < 0) dec = 0;
    } else {
        dec = sig - (int)std::floor(std::log10(val)) - 1; /* sig significant digits */
    }
    return QString::number(val, 'f', dec);
}

PeriodogramTask::PeriodogramTask(lc_compute_ctx_t *ctx, const lc_data_t *data, const lc_periodogram_config_t &cfg, lc_progress_t *progress)
    : m_ctx(ctx), m_data(cloneLcData(data)), m_cfg(cfg), m_progress(progress) {}

PeriodogramTask::~PeriodogramTask() {
    lc_free(&m_data);
}

void PeriodogramTask::run() {
    lc_periodogram_result_t result;
    memset(&result, 0, sizeof(result));
    int rc = lc_compute_periodogram_ctx(m_ctx, &m_data, &m_cfg, &result, m_progress);
    if (rc == 0) {
        QVector<float> nll(result.nll, result.nll + result.nfreq);
        QVector<QPair<double, float>> peaks;
        peaks.reserve(result.npeaks);
        for (int i = 0; i < result.npeaks; ++i)
            peaks.append(qMakePair(result.peaks[i].freq, result.peaks[i].nll));
        double fmin = result.fmin;
        double fstep = result.fstep;
        lc_periodogram_result_free(&result);
        emit finished(fmin, fstep, nll, peaks);
    } else if (rc == 1) {
        emit cancelled();
    } else {
        emit failed(QStringLiteral("Periodogram computation failed"));
    }
}

void PhasedModelWorker::compute(lc_data_t data, lc_periodogram_config_t cfg, double freq) {
    QVector<float> model;
    lc_model_style_t style = LC_MODEL_LINE;
    int rc = -1;
    qint64 ns = 0;
    if (data.n > 0) {
        model.resize((int)data.n);
        QElapsedTimer t;
        t.start();
        rc = lc_compute_phased_model(&data, &cfg, freq, model.data(), &style);
        ns = t.nsecsElapsed();
    }
    lc_free(&data);
    emit done(rc, freq, model, (int)style, ns);
}
