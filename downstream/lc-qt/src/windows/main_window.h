#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QGroupBox>
#include <QVector>
#include <QPair>
#include <QThread>
#include <QTimer>
#include <QElapsedTimer>
#include <functional>

extern "C" {
#include "lc_readout.h"
#include "lc_period.h"
}

#include "customize_period_search.h"

class LightCurvePlotWidget;
class PhasedLightCurveWidget;
class SpectrumPlotWidget;
class ZoomedSpectrumWidget;
class PeriodogramTask;
class PhasedModelWorker;

struct FileEntry {
    QString path;
    int n = 0;
};

struct TargetEntry {
    QString path;
    int cls = 0;
};

class ClassificationDisplay : public QLineEdit {
    Q_OBJECT
    QString *m_labels;
    int m_current;
    QLineEdit *m_typeEdit;
    std::function<void()> m_saveCallback;
public:
    ClassificationDisplay(QString labels[10], QWidget *parent = nullptr);
    void refreshDisplay();
    void setLabels(QString labels[10]);
    void setCurrent(int idx);
    int current() const;
    void setTypeEdit(QLineEdit *typeEdit);
    void setSaveCallback(std::function<void()> cb);

    bool eventFilter(QObject *watched, QEvent *event) override;
};

class PeriodScrollKeyFilter : public QObject {
    Q_OBJECT
    std::function<void(int)> m_start;
    std::function<void()> m_stop;
    std::function<void()> m_zoomIn;
    std::function<void()> m_zoomOut;
public:
    PeriodScrollKeyFilter(std::function<void(int)> start, std::function<void()> stop,
                          std::function<void()> zoomIn, std::function<void()> zoomOut,
                          QObject *parent = nullptr);

    bool eventFilter(QObject *watched, QEvent *event) override;
};

class PeriodModifyKeyFilter : public QObject {
    Q_OBJECT
    std::function<void(double)> m_setFreq;
    std::function<void()> m_toggleDisplay;
    std::function<double()> m_getFreq;

    bool m_multKeyHeld = false;
    bool m_divKeyHeld = false;
    bool m_numberPressed = false;
    double m_baseFreq = 0.0;
public:
    PeriodModifyKeyFilter(std::function<void(double)> setFreq,
                          std::function<void()> toggleDisplay,
                          std::function<double()> getFreq,
                          QObject *parent = nullptr);

    bool eventFilter(QObject *watched, QEvent *event) override;
};

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

    void loadLightCurve(const QString &path);

private:
    void setupUi();
    void setupWiring();
    void loadConfig();
    void saveConfig();

    void launchSpectrum(lc_spec_method_t method);
    void tryDispatchModel();
    void setButtonsEnabled(bool enabled);
    void selectPeakAndAlignExtremum(double peakFreq);
    void moveToPeak(double peakFreq);
    double findPeakNear(double clickFreq, double searchHalfWidth);
    void openPeriodSearchDialog();
    void openPeriodScrollDialog();
    void startStep(int dir);
    void stopStep();

    // Configuration & State
    QString m_labels[10];
    bool m_numpadNav;
    QVector<FileEntry> m_files;
    QVector<TargetEntry> m_targets;
    PeriodSearchSettings m_ps;
    int m_counts[10];
    QString m_currentPath;
    lc_data_t m_lcData;

    // Backend context & progress
    lc_progress_t *m_progress;
    lc_compute_ctx_t *m_computeCtx;

    // Detected peaks & spectrum data
    QVector<QPair<double, float>> m_detectedPeaks;
    double m_specFmin;
    double m_specFstep;
    QVector<float> m_specNll;
    int m_activeMethodIdx;

    // Workers & threads
    QThread *m_workerThread;
    PeriodogramTask *m_task;
    QElapsedTimer m_elapsedTimer;
    QTimer *m_progressTimer;

    QThread *m_modelThread;
    PhasedModelWorker *m_modelWorker;
    bool m_modelBusy;
    bool m_modelRefreshPending;
    double m_pendingModelFreq;

    QTimer *m_stepTimer;
    QElapsedTimer m_stepClock;
    int m_stepDir;

    // UI Widgets
    QLineEdit *m_typeEdit;
    QLineEdit *m_hdj0Edit;
    QLineEdit *m_medEdit;
    QLineEdit *m_madEdit;
    QLineEdit *m_ampEdit;
    QLineEdit *m_objectEdit;
    QLabel *m_noPointsLabel;
    QPushButton *m_displayModelBtn;

    QLineEdit *m_objListEdit;
    QPushButton *m_openObjBtn;
    QLineEdit *m_entryNoEdit;
    QPushButton *m_customizeBtn;
    QPushButton *m_statsBtn;
    ClassificationDisplay *m_classDisplay;

    LightCurvePlotWidget *m_rawPlot;
    PhasedLightCurveWidget *m_phasedPlot;
    SpectrumPlotWidget *m_searchPlot;
    ZoomedSpectrumWidget *m_zoomPlot;

    QPushButton *m_aovBtn;
    QPushButton *m_ihsBtn;
    QPushButton *m_gbBtn;
    QPushButton *m_blsBtn;
    QPushButton *m_stopBtn;
    QPushButton *m_moreBtn;

    QPushButton *m_stepLeftBtn;
    QLineEdit *m_periodVal;
    QPushButton *m_stepRightBtn;
    QPushButton *m_x2Btn;
    QPushButton *m_d2Btn;
    QPushButton *m_otherBtn;
    QMenu *m_otherMenu;

    QPushButton *m_minusBtn;
    QPushButton *m_plusBtn;

    QLabel *m_phasedTimeLabel;
    QLabel *m_progressLabel;

    QAction *m_loadLcAction;
    QAction *m_periodSearchAction;
    QAction *m_periodScrollAction;

    PeriodScrollKeyFilter *m_scrollFilter;
    PeriodModifyKeyFilter *m_modifyFilter;
};

#endif // MAIN_WINDOW_H
