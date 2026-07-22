#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QFrame>
#include <QGroupBox>
#include <QStyle>
#include <QPainter>
#include <QPainterPath>
#include <QKeyEvent>
#include <QtMath>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <functional>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QDir>
#include <QFileDialog>
#include <QMessageBox>
#include <QVector>
#include <QThread>
#include <QTimer>
#include <QStatusBar>
#include <QElapsedTimer>
#include <QAbstractSpinBox>

#include "json.h"

extern "C" {
#include "lc_readout.h"
#include "lc_periodogram.h"
}

#include "windows/customize_labels.h"
#include "windows/classification_stats.h"
#include "windows/lightcurve_plot.h"
#include "windows/phased_lightcurve.h"
#include "windows/spectrum_plot.h"
#include "windows/zoomed_spectrum.h"
#include "windows/customize_period_search.h"

static const char *CONFIG_FILE = ".lc-config.json";

struct FileEntry {
    QString path;
    int n = 0;
};

static void loadConfig(QString labels[10], bool *numpadNav, QVector<FileEntry> *files, PeriodSearchSettings *ps) {
    QFile f(CONFIG_FILE);
    if (!f.open(QIODevice::ReadOnly))
        return;
    QByteArray data = f.readAll();
    f.close();

    json_value_s *root = json_parse(data.constData(), data.size());
    if (!root || root->type != json_type_object) {
        if (root) free(root);
        return;
    }

    json_object_s *obj = json_value_as_object(root);
    for (json_object_element_s *el = obj->start; el; el = el->next) {
        const char *key = el->name->string;

        /* New format: "labels" as JSON array */
        if (strcmp(key, "labels") == 0 && el->value->type == json_type_array) {
            json_array_s *arr = json_value_as_array(el->value);
            int i = 0;
            for (json_array_element_s *ae = arr->start; ae && i < 10; ae = ae->next, i++) {
                if (ae->value->type == json_type_string) {
                    json_string_s *s = json_value_as_string(ae->value);
                    if (s) labels[i] = QString::fromUtf8(s->string, s->string_size);
                }
            }
            continue;
        }

        /* "files" as JSON array of {path, n} objects */
        if (strcmp(key, "files") == 0 && el->value->type == json_type_array) {
            json_array_s *arr = json_value_as_array(el->value);
            for (json_array_element_s *ae = arr->start; ae; ae = ae->next) {
                if (ae->value->type != json_type_object) continue;
                json_object_s *fobj = json_value_as_object(ae->value);
                FileEntry entry;
                for (json_object_element_s *fe = fobj->start; fe; fe = fe->next) {
                    if (strcmp(fe->name->string, "path") == 0 && fe->value->type == json_type_string) {
                        json_string_s *s = json_value_as_string(fe->value);
                        if (s) entry.path = QString::fromUtf8(s->string, s->string_size);
                    } else if (strcmp(fe->name->string, "n") == 0 && fe->value->type == json_type_number) {
                        json_number_s *num = json_value_as_number(fe->value);
                        if (num) entry.n = atoi(num->number);
                    }
                }
                if (!entry.path.isEmpty()) files->append(entry);
            }
            continue;
        }

        /* "period_search" as JSON object */
        if (strcmp(key, "period_search") == 0 && el->value->type == json_type_object && ps) {
            json_object_s *pobj = json_value_as_object(el->value);
            for (json_object_element_s *pe = pobj->start; pe; pe = pe->next) {
                if (pe->value->type != json_type_number) continue;
                json_number_s *num = json_value_as_number(pe->value);
                if (!num) continue;
                const char *pk = pe->name->string;
                if (strcmp(pk, "nterms") == 0) ps->nterms = atoi(num->number);
                else if (strcmp(pk, "oversampling") == 0) ps->oversampling = atof(num->number);
                else if (strcmp(pk, "fmin") == 0) ps->fmin = atof(num->number);
                else if (strcmp(pk, "fmax") == 0) ps->fmax = atof(num->number);
                else if (strcmp(pk, "zoom_factor") == 0) ps->zoomFactor = atof(num->number);
                else if (strcmp(pk, "search_radius") == 0) ps->searchRadius = atof(num->number);
                else if (strcmp(pk, "pswf") == 0) ps->pswf = atoi(num->number);
                else if (strcmp(pk, "oversmoothing") == 0) ps->oversmoothing = atof(num->number);
                else if (strcmp(pk, "nbins") == 0) ps->nbins = atoi(num->number);
            }
            continue;
        }

        if (el->value->type == json_type_string) {
            json_string_s *s = json_value_as_string(el->value);
            if (!s) continue;
            QString val = QString::fromUtf8(s->string, s->string_size);

            if (strcmp(key, "numpad_nav") == 0) {
                *numpadNav = (val == "true");
            }

            /* Legacy format: "label0"..."label9" */
            for (int i = 0; i < 10; ++i) {
                char lkey[8];
                snprintf(lkey, sizeof(lkey), "label%d", i);
                if (strcmp(key, lkey) == 0) {
                    labels[i] = val;
                    break;
                }
            }
        }
    }
    free(root);
}

static QString escapeJsonString(const QString &s) {
    QString out = s;
    out.replace('\\', "\\\\").replace('"', "\\\"");
    return out;
}

static void saveConfig(const QString labels[10], bool numpadNav, const QVector<FileEntry> &files, const PeriodSearchSettings &ps) {
    QString json = "{\"labels\":[";
    for (int i = 0; i < 10; ++i) {
        json += QString("\"%1\"").arg(escapeJsonString(labels[i]));
        if (i < 9) json += ',';
    }
    json += QString("],\"numpad_nav\":\"%1\"").arg(numpadNav ? "true" : "false");

    json += ",\"files\":[";
    for (int i = 0; i < files.size(); ++i) {
        json += QString("{\"path\":\"%1\",\"n\":%2}").arg(escapeJsonString(files[i].path)).arg(files[i].n);
        if (i < files.size() - 1) json += ',';
    }
    json += "]";

    json += QString(",\"period_search\":{"
                    "\"nterms\":%1,"
                    "\"oversampling\":%2,"
                    "\"fmin\":%3,"
                    "\"fmax\":%4,"
                    "\"zoom_factor\":%5,"
                    "\"search_radius\":%6,"
                    "\"pswf\":%7,"
                    "\"oversmoothing\":%8,"
                    "\"nbins\":%9}")
                .arg(ps.nterms)
                .arg(ps.oversampling, 0, 'g', 10)
                .arg(ps.fmin, 0, 'g', 10)
                .arg(ps.fmax, 0, 'g', 10)
                .arg(ps.zoomFactor, 0, 'g', 10)
                .arg(ps.searchRadius, 0, 'g', 10)
                .arg(ps.pswf)
                .arg(ps.oversmoothing, 0, 'g', 10)
                .arg(ps.nbins);
    json += "}";

    QFile f(CONFIG_FILE);
    if (!f.open(QIODevice::WriteOnly))
        return;
    f.write(json.toUtf8());
    f.close();
}

class MockPlotWidget : public QFrame {
    QString m_title;
public:
    MockPlotWidget(const QString &title, bool isPhased = false, bool isSearch = false, bool isModify = false, QWidget *parent = nullptr)
        : QFrame(parent), m_title(title) {
        setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
        setMinimumHeight(150);
        
        // Graphite dark plot background
        setStyleSheet("background-color: #1a1a1f; border: 1px solid #32323d; border-radius: 4px;");
    }

protected:
    void paintEvent(QPaintEvent *event) override {
        QFrame::paintEvent(event);
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing);

        int w = width();
        int h = height();
        int padding = 35;

        // Draw soft dark grid
        painter.setPen(QPen(QColor(42, 42, 53), 1, Qt::DashLine));
        int gridCount = 6;
        for (int i = 1; i < gridCount; ++i) {
            // Horizontal lines
            int y = padding + i * (h - 2 * padding) / gridCount;
            painter.drawLine(padding, y, w - padding, y);
            // Vertical lines
            int x = padding + i * (w - 2 * padding) / gridCount;
            painter.drawLine(x, padding, x, h - padding);
        }

        // Draw clean axes
        painter.setPen(QPen(QColor(100, 116, 139), 1.2));
        painter.drawLine(padding, h - padding, w - padding, h - padding); // X-axis
        painter.drawLine(padding, padding, padding, h - padding);       // Y-axis

        // Draw title
        painter.setPen(QColor(203, 213, 225));
        painter.setFont(QFont("sans-serif", 9, QFont::Bold));
        painter.drawText(QRect(padding, 4, w - 2 * padding, 20), Qt::AlignCenter, m_title);
    }
};

class ClassificationDisplay : public QLineEdit {
    QString *m_labels;
    int m_current;
    QLineEdit *m_typeEdit;
public:
    ClassificationDisplay(QString labels[10], QWidget *parent = nullptr)
        : QLineEdit(parent), m_labels(labels), m_current(0), m_typeEdit(nullptr) {
        setReadOnly(true);
        setFocusPolicy(Qt::NoFocus);
        setMinimumWidth(200);
        setAlignment(Qt::AlignCenter);
        refreshDisplay();
    }

    void refreshDisplay() {
        setText(QString("%1 — %2").arg(m_current).arg(m_labels[m_current]));
    }

    void setLabels(QString labels[10]) {
        m_labels = labels;
        refreshDisplay();
    }

    void setCurrent(int idx) {
        m_current = idx;
        refreshDisplay();
    }

    void setTypeEdit(QLineEdit *typeEdit) {
        m_typeEdit = typeEdit;
    }

    bool eventFilter(QObject *watched, QEvent *event) override {
        if (event->type() == QEvent::KeyPress) {
            QKeyEvent *ke = static_cast<QKeyEvent *>(event);
            int key = ke->key();
            if (key >= Qt::Key_0 && key <= Qt::Key_9) {
                QWidget *focus = QApplication::focusWidget();
                bool textHasFocus = focus && (qobject_cast<QLineEdit *>(focus) || qobject_cast<QAbstractSpinBox *>(focus));
                if (!textHasFocus) {
                    m_current = key - Qt::Key_0;
                    refreshDisplay();
                    return true;
                }
            } else if (key == Qt::Key_Return || key == Qt::Key_Enter) {
                QWidget *focus = QApplication::focusWidget();
                bool textHasFocus = focus && (qobject_cast<QLineEdit *>(focus) || qobject_cast<QAbstractSpinBox *>(focus));
                if (!textHasFocus) {
                    if (m_typeEdit) {
                        m_typeEdit->setText(m_labels[m_current]);
                    }
                    return true;
                }
            }
        }
        return QObject::eventFilter(watched, event);
    }
};


/* Deep copy of an lc_data_t for use on the worker thread (frees only x/y/dy; _buf stays null). */
static lc_data_t cloneLcData(const lc_data_t *src) {
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

/* Worker that runs the (blocking, multithreaded) periodogram computation off the GUI thread. */
class PeriodogramTask : public QObject {
    Q_OBJECT
public:
    PeriodogramTask(lc_compute_ctx_t *ctx, const lc_data_t *data, const lc_periodogram_config_t &cfg, lc_progress_t *progress)
        : m_ctx(ctx), m_data(cloneLcData(data)), m_cfg(cfg), m_progress(progress) {}
    ~PeriodogramTask() override { lc_free(&m_data); }

public slots:
    void run() {
        lc_periodogram_result_t result;
        memset(&result, 0, sizeof(result));
        int rc = lc_compute_periodogram_ctx(m_ctx, &m_data, &m_cfg, &result, m_progress);
        if (rc == 0) {
            QVector<float> nll(result.nll, result.nll + result.nfreq);
            double fmin = result.fmin;
            double fstep = result.fstep;
            lc_periodogram_result_free(&result);
            emit finished(fmin, fstep, nll);
        } else if (rc == 1) {
            emit cancelled();
        } else {
            emit failed(QStringLiteral("Periodogram computation failed"));
        }
    }

signals:
    void finished(double fmin, double fstep, QVector<float> nll);
    void cancelled();
    void failed(QString msg);

private:
    lc_compute_ctx_t *m_ctx;
    lc_data_t m_data;
    lc_periodogram_config_t m_cfg;
    lc_progress_t *m_progress;
};

/* Format the period (1/freq) the same way upstream ihsnpeaks formats peak
   periods with --period (src/utils/strings.h: custom_coordinate_dtoa):
   significant digits = precision + ceil(log10(freq)) clamped to [1, 17], where
   precision = 1 + log10(DeltaT * oversampling * nterms); the fractional part is
   zero-padded so the value shows that many significant digits. */
static QString formatPeriodUpstream(double freq, double deltaT, double oversampling, int nterms) {
    if (!(freq > 0.0) || !std::isfinite(freq)) return QStringLiteral("inf");
    double arg = deltaT * oversampling * (double)nterms;
    int precision = (arg > 0.0) ? 1 + (int)std::log10(arg) : 1;
    int sig = precision + (int)std::ceil(std::log10(freq));
    if (sig < 1) sig = 1;
    if (sig > 17) sig = 17;
    double period = 1.0 / freq;
    int dec;
    if (period >= 1.0) {
        int intDigits = (int)std::ceil(std::log10(period));
        if (intDigits < 0) intDigits = 0;
        dec = sig - intDigits;
        if (dec < 0) dec = 0;
    } else {
        dec = sig - (int)std::floor(std::log10(period)) - 1; /* sig significant digits */
    }
    return QString::number(period, 'f', dec);
}


int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    
    // Set classic Motif-like color palette and style hint
    app.setStyle("Fusion");
    
    QMainWindow window;
    window.setWindowTitle("lc-qt v ?.?? 21/07/2026");
    window.resize(900, 850);

    // Classification state
    QString labels[10] = {"nonvar", "var", "unknown", "unknown", "unknown",
                          "unknown", "unknown", "unknown", "unknown", "unknown"};
    bool numpadNav = false;
    QVector<FileEntry> files;
    PeriodSearchSettings ps;
    loadConfig(labels, &numpadNav, &files, &ps);
    int counts[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    // Menu Bar
    QMenuBar *menuBar = window.menuBar();
    QMenu *fileMenu = menuBar->addMenu("File");
    QAction *loadLcAction = fileMenu->addAction("Load Light Curve...");
    fileMenu->addAction("Save Light Curve");
    fileMenu->addAction("Save Light Curve As...");
    fileMenu->addAction("Open Object List...");
    fileMenu->addAction("Save Object List...");
    fileMenu->addSeparator();
    fileMenu->addAction("Quit", &app, &QCoreApplication::quit);
    
    QMenu *optionsMenu = menuBar->addMenu("Options");
    //optionsMenu->addAction("Customize In/Out...");
    //optionsMenu->addAction("Preview");
    //optionsMenu->addAction("Classes");
    //optionsMenu->addAction("Frequency analysis");
    optionsMenu->addAction("Customize Period Search...");
    QAction *periodSearchAction = optionsMenu->actions().last();
    //optionsMenu->addAction("Light Curve");
    optionsMenu->addAction("Plot Options");
    optionsMenu->addAction("Period Scroll");
    //optionsMenu->addAction("Multi-Period Window");
    //optionsMenu->addAction("Auto Equalize");
    optionsMenu->addAction("Auto Power Spec Calc");
    //optionsMenu->addAction("Auto Adjust Spectrum Max");
    
    menuBar->addAction("Help");
    
    // Main layout
    QWidget *centralWidget = new QWidget(&window);
    window.setCentralWidget(centralWidget);
    QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);
    mainLayout->setContentsMargins(6, 6, 6, 6);
    mainLayout->setSpacing(6);
    
    // Global style sheet for a comfortable graphite dark palette
    centralWidget->setStyleSheet(
        "QWidget { background-color: #242429; color: #e2e8f0; font-family: -apple-system, BlinkMacSystemFont, Roboto, 'Segoe UI', sans-serif; font-size: 11px; }"
        "QMainWindow { background-color: #242429; }"
        "QMenuBar { background-color: #1c1c20; color: #cbd5e1; border-bottom: 1px solid #2d2d35; }"
        "QMenuBar::item:selected { background-color: #3b82f6; color: #ffffff; }"
        "QMenu { background-color: #1c1c20; color: #cbd5e1; border: 1px solid #2d2d35; }"
        "QMenu::item:selected { background-color: #3b82f6; color: #ffffff; }"
        "QLineEdit { background-color: #17171b; border: 1px solid #32323b; border-radius: 4px; color: #ffffff; padding: 3px; }"
        "QLineEdit:focus { border: 1px solid #3b82f6; }"
        "QPushButton { background-color: #303038; border: 1px solid #3d3d47; border-radius: 4px; color: #f1f5f9; padding: 4px 8px; font-weight: bold; }"
        "QPushButton:hover { background-color: #3b82f6; border-color: #60a5fa; }"
        "QPushButton:pressed { background-color: #1d4ed8; }"
        "QGroupBox { border: 1px solid #32323b; border-radius: 6px; margin-top: 10px; padding: 6px; font-weight: bold; color: #60a5fa; }"
        "QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; background-color: #242429; }"
        "QLabel { color: #cbd5e1; }"
    );



    
    // Type / HJDO / Max,Amp / Object bar
    QHBoxLayout *topBarLayout = new QHBoxLayout();
    
    topBarLayout->addWidget(new QLabel("Type"));
    QLineEdit *typeEdit = new QLineEdit("unknown");
    typeEdit->setMaximumWidth(120);
    typeEdit->setReadOnly(true);
    topBarLayout->addWidget(typeEdit);
    
    topBarLayout->addWidget(new QLabel("HJDO"));
    QLineEdit *hjdoEdit = new QLineEdit("2450000");
    hjdoEdit->setMaximumWidth(80);
    topBarLayout->addWidget(hjdoEdit);
    
    topBarLayout->addWidget(new QLabel("Max"));
    QLineEdit *maxEdit = new QLineEdit("0");
    maxEdit->setMaximumWidth(60);
    topBarLayout->addWidget(maxEdit);

    topBarLayout->addWidget(new QLabel("Amp"));
    QLineEdit *ampEdit = new QLineEdit("0");
    ampEdit->setMaximumWidth(60);
    topBarLayout->addWidget(ampEdit);
    
    topBarLayout->addWidget(new QLabel("Object"));
    QLineEdit *objectEdit = new QLineEdit();
    topBarLayout->addWidget(objectEdit);
    
    mainLayout->addLayout(topBarLayout);
    
    // Raw and Phased Plots Row
    QHBoxLayout *plotsLayout = new QHBoxLayout();
    LightCurvePlotWidget *rawPlot = new LightCurvePlotWidget("Raw Light Curve");
    PhasedLightCurveWidget *phasedPlot = new PhasedLightCurveWidget("Phased Light Curve");
    plotsLayout->addWidget(rawPlot);
    plotsLayout->addWidget(phasedPlot);
    mainLayout->addLayout(plotsLayout, 3); // stretch factor 3

    
    // Object List / Light Curve row
    QHBoxLayout *listAndCurveLayout = new QHBoxLayout();
    
    QGroupBox *objGroupBox = new QGroupBox("Object List");
    QHBoxLayout *objLayout = new QHBoxLayout(objGroupBox);
    QLineEdit *objListEdit = new QLineEdit();
    QPushButton *openObjBtn = new QPushButton("Open");
    QLabel *entryNoLabel = new QLabel("Entry No.");
    QLineEdit *entryNoEdit = new QLineEdit("0");
    entryNoEdit->setMaximumWidth(50);
    objLayout->addWidget(objListEdit);
    objLayout->addWidget(openObjBtn);
    objLayout->addWidget(entryNoLabel);
    objLayout->addWidget(entryNoEdit);
    listAndCurveLayout->addWidget(objGroupBox, 1);
    
    QGroupBox *lcGroupBox = new QGroupBox("Light Curve");
    QHBoxLayout *lcLayout = new QHBoxLayout(lcGroupBox);
    QLabel *noPointsLabel = new QLabel("No. points    0");
    QPushButton *equBtn = new QPushButton("Equ");
    lcLayout->addWidget(noPointsLabel);
    lcLayout->addStretch();
    lcLayout->addWidget(equBtn);
    listAndCurveLayout->addWidget(lcGroupBox, 1);
    
    mainLayout->addLayout(listAndCurveLayout);
    
    // Log Files Section — split into Statistics (left) + Classification (right)
    QHBoxLayout *logFilesLayout = new QHBoxLayout();

    QGroupBox *statsGroupBox = new QGroupBox("Statistics");
    QHBoxLayout *statsLayout = new QHBoxLayout(statsGroupBox);
    QPushButton *statsBtn = new QPushButton("Class stats");
    statsLayout->addStretch();
    statsLayout->addWidget(statsBtn);
    statsLayout->addStretch();
    logFilesLayout->addWidget(statsGroupBox, 1);

    QGroupBox *classGroupBox = new QGroupBox("Classification");
    QHBoxLayout *classLayout = new QHBoxLayout(classGroupBox);

    QPushButton *customizeBtn = new QPushButton("Customize labels");
    ClassificationDisplay *classDisplay = new ClassificationDisplay(labels);
    classDisplay->setTypeEdit(typeEdit);

    classLayout->addStretch();
    classLayout->addWidget(customizeBtn);
    classLayout->addStretch();
    classLayout->addWidget(new QLabel("Current:"));
    classLayout->addWidget(classDisplay);
    classLayout->addStretch();
    logFilesLayout->addWidget(classGroupBox, 1);

    mainLayout->addLayout(logFilesLayout);

    app.installEventFilter(classDisplay);

    // Light curve loading logic
    std::function<void()> onLightCurveLoaded;

    lc_data_t lcData;
    memset(&lcData, 0, sizeof(lcData));

    auto loadLightCurve = [&](const QString &path) {
        if (path.isEmpty()) return;

        QFileInfo fi(path);
        QString suffix = fi.suffix().toLower();

        lc_data_t newData;
        int rc;
        if (suffix == "fits") {
            rc = lc_load_fits(path.toUtf8().constData(), &newData);
        } else {
            rc = lc_load_dat(path.toUtf8().constData(), &newData);
        }

        if (rc != 0 || newData.n == 0) {
            QMessageBox::warning(&window, "Load Error",
                                 QString("Failed to load light curve from:\n%1").arg(path));
            return;
        }

        /* Detrend: remove linear trend, add mean back */
        lc_detrend(&newData);

        /* Free previous data and take ownership of new */
        lc_free(&lcData);
        lcData = newData;

        /* Update plot */
        rawPlot->setData(&lcData);
        phasedPlot->setData(&lcData);

        /* Update UI labels */
        noPointsLabel->setText(QString("No. points    %1").arg(lcData.n));
        objectEdit->setText(fi.fileName());

        /* Update files list in config */
        bool found = false;
        for (int i = 0; i < files.size(); i++) {
            if (files[i].path == path) {
                files[i].n = (int)lcData.n;
                found = true;
                break;
            }
        }
        if (!found) {
            files.append(FileEntry{path, (int)lcData.n});
        }
        saveConfig(labels, numpadNav, files, ps);
        if (onLightCurveLoaded) onLightCurveLoaded();
    };

    QObject::connect(loadLcAction, &QAction::triggered, [&]() {
        QString path = QFileDialog::getOpenFileName(
            &window, "Load Light Curve", QString(),
            "Light Curve Files (*.dat *.fits);;All Files (*)");
        loadLightCurve(path);
    });

    QObject::connect(customizeBtn, &QPushButton::clicked, [&]() {
        bool reopen = true;
        while (reopen) {
            reopen = false;
            CustomizeLabelsDialog dlg(labels, &numpadNav, &window);
            QObject::connect(&dlg, &CustomizeLabelsDialog::labelsChanged, classDisplay, [classDisplay, labels, &numpadNav, &files, &ps]() {
                classDisplay->refreshDisplay();
                saveConfig(labels, numpadNav, files, ps);
            });
            dlg.exec();
            if (dlg.shouldReopen()) {
                reopen = true;
            }
        }
    });

    QObject::connect(statsBtn, &QPushButton::clicked, [&]() {
        ClassificationStatsDialog dlg(labels, counts, &window);
        dlg.exec();
    });

    // Period Search / Period Modify section
    QHBoxLayout *periodLayout = new QHBoxLayout();
    
    QGroupBox *searchGroupBox = new QGroupBox("Period Search");
    QVBoxLayout *searchLayout = new QVBoxLayout(searchGroupBox);
    QHBoxLayout *searchBtns = new QHBoxLayout();
    QPushButton *aovBtn = new QPushButton("AoV");
    QPushButton *ihsBtn = new QPushButton("IHS");
    QPushButton *gbBtn = new QPushButton("GB");
    QPushButton *blsBtn = new QPushButton("BLS");
    QPushButton *stopBtn = new QPushButton("Stop");
    QPushButton *moreBtn = new QPushButton("…");
    searchBtns->addWidget(aovBtn);
    searchBtns->addWidget(ihsBtn);
    searchBtns->addWidget(gbBtn);
    searchBtns->addWidget(blsBtn);
    searchBtns->addWidget(stopBtn);
    searchBtns->addWidget(moreBtn);
    searchLayout->addLayout(searchBtns);
    SpectrumPlotWidget *searchPlot = new SpectrumPlotWidget("Negative Log-Likelihood");
    searchLayout->addWidget(searchPlot, 1);
    periodLayout->addWidget(searchGroupBox, 1);

    
    QGroupBox *modifyGroupBox = new QGroupBox("Period Modify");
    QVBoxLayout *modifyLayout = new QVBoxLayout(modifyGroupBox);
    QHBoxLayout *modifyBtns = new QHBoxLayout();
    QPushButton *stepLeftBtn = new QPushButton("<");
    QLineEdit *periodVal = new QLineEdit("0.000");
    periodVal->setMaximumWidth(80);
    QPushButton *stepRightBtn = new QPushButton(">");
    QPushButton *x2Btn = new QPushButton("x2");
    QPushButton *d2Btn = new QPushButton("/2");
    QPushButton *otherBtn = new QPushButton("Other");
    QMenu *otherMenu = new QMenu(otherBtn);
    otherBtn->setMenu(otherMenu);
    modifyBtns->addWidget(stepLeftBtn);
    modifyBtns->addWidget(periodVal);
    modifyBtns->addWidget(stepRightBtn);
    modifyBtns->addWidget(x2Btn);
    modifyBtns->addWidget(d2Btn);
    modifyBtns->addWidget(otherBtn);
    modifyBtns->addWidget(new QPushButton("◄")); /* mock-up: wired once all four widgets work */
    modifyBtns->addWidget(new QPushButton("►")); /* mock-up: wired once all four widgets work */
    modifyLayout->addLayout(modifyBtns);
    
    // Bottom layout for zoomed spectrum showing [-] and [+] on the side
    QHBoxLayout *bottomModifyArea = new QHBoxLayout();
    QPushButton *minusBtn = new QPushButton("-");
    minusBtn->setMaximumSize(25, 50);
    QPushButton *plusBtn = new QPushButton("+");
    plusBtn->setMaximumSize(25, 50);
    
    ZoomedSpectrumWidget *zoomPlot = new ZoomedSpectrumWidget("Zoomed Spectrum");

    
    QVBoxLayout *leftCtrl = new QVBoxLayout();
    leftCtrl->addWidget(minusBtn);
    leftCtrl->addStretch();
    
    QVBoxLayout *rightCtrl = new QVBoxLayout();
    rightCtrl->addWidget(plusBtn);
    rightCtrl->addStretch();
    
    bottomModifyArea->addLayout(leftCtrl);
    bottomModifyArea->addWidget(zoomPlot, 1);
    bottomModifyArea->addLayout(rightCtrl);
    
    modifyLayout->addLayout(bottomModifyArea, 1);
    
    periodLayout->addWidget(modifyGroupBox, 1);
    
    mainLayout->addLayout(periodLayout, 2); // stretch factor 2

    // ---- Zoomed spectrum wiring ----
    QObject::connect(minusBtn, &QPushButton::clicked, zoomPlot, &ZoomedSpectrumWidget::zoomOut);
    QObject::connect(plusBtn, &QPushButton::clicked, zoomPlot, &ZoomedSpectrumWidget::zoomIn);
    QObject::connect(searchPlot, &SpectrumPlotWidget::frequencyClicked, zoomPlot,
                     [zoomPlot](double freq) { zoomPlot->selectFrequency(freq, true); });
    QObject::connect(zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged,
                     searchPlot, &SpectrumPlotWidget::setSelectedFrequency);
    QObject::connect(zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged,
                     phasedPlot, &PhasedLightCurveWidget::setFrequency);
    QObject::connect(searchPlot, &SpectrumPlotWidget::rangeSelected,
                     zoomPlot, &ZoomedSpectrumWidget::setViewFromSelection);

    /* x2 / /2 buttons and the "Other" multiplier menu (x3 /3 x5 /5 x7 /7):
       scale the selected (pivot) frequency; the FOV re-centres on the pivot. */
    static const struct {
        const char *label;
        double factor;
    } otherMults[] = {{"x3", 3.0}, {"/3", 1.0 / 3.0}, {"x5", 5.0}, {"/5", 1.0 / 5.0}, {"x7", 7.0}, {"/7", 1.0 / 7.0}};
    for (const auto &m : otherMults) {
        QAction *act = otherMenu->addAction(QString::fromLatin1(m.label));
        QObject::connect(act, &QAction::triggered, [zoomPlot, m] { zoomPlot->multiplyFrequency(m.factor); });
    }
    QObject::connect(x2Btn, &QPushButton::clicked, [zoomPlot] { zoomPlot->multiplyFrequency(2.0); });
    QObject::connect(d2Btn, &QPushButton::clicked, [zoomPlot] { zoomPlot->multiplyFrequency(0.5); });

    /* Period field: auto-displays the period (1/f) of the pivot formatted like
       upstream ihsnpeaks --period peak output; editing it moves the pivot. */
    QObject::connect(zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged, periodVal,
                     [&lcData, &ps, periodVal](double freq) {
                         double dT = (lcData.n > 1) ? lcData.x[lcData.n - 1] - lcData.x[0] : 0.0;
                         periodVal->setText(formatPeriodUpstream(freq, dT, ps.oversampling, ps.nterms));
                     });
    QObject::connect(periodVal, &QLineEdit::editingFinished, [periodVal, zoomPlot]() {
        bool ok = false;
        double p = periodVal->text().toDouble(&ok);
        if (ok && p > 0.0) zoomPlot->selectFrequency(1.0 / p, false);
    });

    /* Arrow buttons: while held, move the selected frequency at 1/DeltaT per
       second (DeltaT = time from the first to the last measurement). */
    QTimer *stepTimer = new QTimer(&window);
    stepTimer->setInterval(50);
    QElapsedTimer stepClock;
    int stepDir = 0;
    QObject::connect(stepTimer, &QTimer::timeout, [&]() {
        double dt = (double)stepClock.restart() * 1e-9; /* ns -> s */
        double dT = (lcData.n > 1) ? lcData.x[lcData.n - 1] - lcData.x[0] : 0.0;
        if (dT > 0.0 && stepDir != 0) zoomPlot->stepFrequency((double)stepDir * dt / dT);
    });
    auto startStep = [&](int dir) {
        stepDir = dir;
        stepClock.start();
        stepTimer->start();
    };
    auto stopStep = [&]() {
        stepTimer->stop();
        stepDir = 0;
    };
    QObject::connect(stepLeftBtn, &QPushButton::pressed, [&] { startStep(-1); });
    QObject::connect(stepLeftBtn, &QPushButton::released, stopStep);
    QObject::connect(stepRightBtn, &QPushButton::pressed, [&] { startStep(+1); });
    QObject::connect(stepRightBtn, &QPushButton::released, stopStep);

    // ---- Periodogram spectrum computation wiring ----
    lc_progress_t *progress = lc_progress_create();
    lc_compute_ctx_t *computeCtx = lc_compute_ctx_create(0 /* auto: physical cores */);

    QLabel *progressLabel = new QLabel();
    progressLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    window.statusBar()->addPermanentWidget(progressLabel);

    QThread *workerThread = nullptr;
    PeriodogramTask *task = nullptr;
    QElapsedTimer elapsedTimer;
    QTimer *progressTimer = new QTimer(&window);
    progressTimer->setInterval(100);

    auto setButtonsEnabled = [&](bool enabled) {
        bool hasData = (lcData.n > 0);
        aovBtn->setEnabled(enabled && hasData);
        ihsBtn->setEnabled(enabled && hasData);
        gbBtn->setEnabled(enabled && hasData);
        blsBtn->setEnabled(enabled && hasData);
        stopBtn->setEnabled(!enabled);
    };
    onLightCurveLoaded = [&]() { setButtonsEnabled(true); };

    QObject::connect(progressTimer, &QTimer::timeout, [&]() {
        uint32_t done = lc_progress_done(progress);
        uint32_t total = lc_progress_total(progress);
        if (total == 0) {
            searchPlot->setProgressText("...");
            progressLabel->setText("Computation in progress...");
            return;
        }
        double pct = 100.0 * (double)done / (double)total;
        searchPlot->setProgressText(QString("%1%").arg(pct, 0, 'f', 1));
        QString timeLeft;
        double elapsed = (double)elapsedTimer.elapsed() / 1000.0;
        if (elapsed > 0.001 && done > 0) {
            double remaining = (double)(total - done) * elapsed / (double)done;
            if (remaining >= 60.0) {
                int mins = (int)(remaining / 60.0);
                int secs = (int)(remaining - (double)mins * 60.0);
                timeLeft = QString("%1m %2s left").arg(mins).arg(secs);
            } else {
                timeLeft = QString("%1s left").arg(remaining, 0, 'f', 2);
            }
        } else {
            timeLeft = "calculating...";
        }
        progressLabel->setText(QString("Computation in progress: %1% complete | %2").arg(pct, 0, 'f', 1).arg(timeLeft));
    });

    auto launchSpectrum = [&](lc_spec_method_t method) {
        if (lcData.n == 0) {
            QMessageBox::warning(&window, "No Data", "Load a light curve first.");
            return;
        }
        if (workerThread) return; /* already running */

        lc_periodogram_config_t cfg;
        memset(&cfg, 0, sizeof(cfg));
        cfg.method = method;
        cfg.nterms = ps.nterms;
        cfg.oversampling = ps.oversampling;
        cfg.fmin = ps.fmin;
        cfg.fmax = ps.fmax;
        cfg.pswf = ps.pswf;
        cfg.nthreads = 0; /* auto: all online CPUs */
        cfg.oversmoothing = ps.oversmoothing;
        cfg.nbins = ps.nbins;

        lc_progress_reset(progress);
        setButtonsEnabled(false);
        elapsedTimer.start();
        searchPlot->setProgressText("...");
        progressLabel->setText("Computation in progress...");
        progressTimer->start();

        workerThread = new QThread();
        task = new PeriodogramTask(computeCtx, &lcData, cfg, progress);
        task->moveToThread(workerThread);

        QObject::connect(workerThread, &QThread::started, task, &PeriodogramTask::run);
        QObject::connect(task, &PeriodogramTask::finished, searchPlot,
                         [searchPlot](double fmin, double fstep, QVector<float> nll) { searchPlot->setData(fmin, fstep, nll); });
        QObject::connect(task, &PeriodogramTask::finished, zoomPlot,
                         [zoomPlot, &ps](double fmin, double fstep, QVector<float> nll) {
                             zoomPlot->setZoomFactor(ps.zoomFactor);
                             zoomPlot->setFullSpectrum(fmin, fstep, nll);
                         });
        QObject::connect(task, &PeriodogramTask::failed, &window,
                         [&window](QString msg) { QMessageBox::warning(&window, "Computation Failed", msg); });

        auto onComplete = [&]() {
            progressTimer->stop();
            progressLabel->clear();
            searchPlot->setProgressText("");
            setButtonsEnabled(true);
            if (workerThread) workerThread->quit();
        };
        QObject::connect(task, &PeriodogramTask::finished, &window, [onComplete](double, double, QVector<float>) { onComplete(); });
        QObject::connect(task, &PeriodogramTask::cancelled, &window, [onComplete]() { onComplete(); });
        QObject::connect(task, &PeriodogramTask::failed, &window, [onComplete](QString) { onComplete(); });

        QObject::connect(workerThread, &QThread::finished, task, &QObject::deleteLater);
        QObject::connect(workerThread, &QThread::finished, workerThread, &QObject::deleteLater);
        QObject::connect(workerThread, &QThread::finished, &window, [&]() {
            workerThread = nullptr;
            task = nullptr;
        });

        workerThread->start();
    };

    QObject::connect(aovBtn, &QPushButton::clicked, [&]() { launchSpectrum(LC_SPEC_AOV); });
    QObject::connect(ihsBtn, &QPushButton::clicked, [&]() { launchSpectrum(LC_SPEC_IHS); });
    QObject::connect(gbBtn, &QPushButton::clicked, [&]() { launchSpectrum(LC_SPEC_GB); });
    QObject::connect(blsBtn, &QPushButton::clicked, [&]() { launchSpectrum(LC_SPEC_BLS); });
    QObject::connect(stopBtn, &QPushButton::clicked, [&]() { lc_progress_request_cancel(progress); });

    auto openPeriodSearchDialog = [&]() {
        CustomizePeriodSearchDialog dlg(&ps, &window);
        if (dlg.exec() == QDialog::Accepted) {
            saveConfig(labels, numpadNav, files, ps);
        }
    };
    QObject::connect(moreBtn, &QPushButton::clicked, openPeriodSearchDialog);
    QObject::connect(periodSearchAction, &QAction::triggered, openPeriodSearchDialog);

    setButtonsEnabled(true); /* initial state: disabled until a light curve is loaded */

    window.show();

    /* Auto-load file passed as argv[1] */
    if (argc > 1) {
        loadLightCurve(QString::fromLocal8Bit(argv[1]));
    }

    int rc = app.exec();
    lc_compute_ctx_destroy(computeCtx);
    lc_progress_destroy(progress);
    return rc;
}

#include "main.moc"
