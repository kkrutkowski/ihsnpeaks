#pragma once

#include <QDialog>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QLineEdit>
#include <QCheckBox>

/* Persisted period-search configuration (mirrors the ihsnpeaks CLI options). */
struct PeriodSearchSettings {
    int nterms = 3;            /* harmonics (-d) */
    double oversampling = 5.0; /* oversampling factor (-o) */
    double fmin = 0.0;         /* min frequency (-f); 0 = auto */
    double fmax = 24.0;        /* max frequency (positional) */
    double zoomFactor = 4.0;   /* zoomed spectrum: initial zoom on full-spectrum click */
    double searchRadius = 0.1; /* stored; no ihsnpeaks spectrum mapping */
    int pswf = 43;             /* NuFFT backend: 43 (pswf43) or 21 (pswf21) */
    double oversmoothing = 0.2;/* -> gbAlpha (GBLS) AND blsMinRelWidth (BLS) */
    int nbins = 10;            /* -> blsWidthCount (BLS only) */
    bool autoCenter = true;    /* automatically center on highest peak after computation */
};

class CustomizePeriodSearchDialog : public QDialog {
    Q_OBJECT
public:
    explicit CustomizePeriodSearchDialog(PeriodSearchSettings *settings, QWidget *parent = nullptr);

private slots:
    void apply();

private:
    PeriodSearchSettings *m_settings;
    QSpinBox *m_terms;
    QDoubleSpinBox *m_oversampling;
    QLineEdit *m_fmin;
    QDoubleSpinBox *m_fmax;
    QDoubleSpinBox *m_zoom;
    QDoubleSpinBox *m_radius;
    QComboBox *m_pswf;
    QDoubleSpinBox *m_oversmoothing;
    QSpinBox *m_nbins;
    QCheckBox *m_autoCenter;
};
