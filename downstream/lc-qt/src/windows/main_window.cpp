#include "main_window.h"

#include <QApplication>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QFile>
#include <QFileInfo>
#include <QFileDialog>
#include <QMessageBox>
#include <QStatusBar>
#include <QKeyEvent>
#include <QAbstractSpinBox>
#include <cmath>
#include <array>
#include <algorithm>

#include "json.h"
#include "period_worker.h"
#include "customize_labels.h"
#include "classification_stats.h"
#include "lightcurve_plot.h"
#include "phased_lightcurve.h"
#include "spectrum_plot.h"
#include "zoomed_spectrum.h"
#include "period_scroll.h"

static const char *CONFIG_FILE = ".lc-config.json";

static QString escapeJsonString(const QString &s) {
    QString out = s;
    out.replace('\\', "\\\\").replace('"', "\\\"");
    return out;
}

// ---------------------------------------------------------------------------
// ClassificationDisplay
// ---------------------------------------------------------------------------

ClassificationDisplay::ClassificationDisplay(QString labels[10], QWidget *parent)
    : QLineEdit(parent), m_labels(labels), m_current(0), m_typeEdit(nullptr) {
    setReadOnly(true);
    setFocusPolicy(Qt::NoFocus);
    setMinimumWidth(120);
    setAlignment(Qt::AlignCenter);
    refreshDisplay();
}

void ClassificationDisplay::refreshDisplay() {
    setText(QString("%1 — %2").arg(m_current).arg(m_labels[m_current]));
}

void ClassificationDisplay::setLabels(QString labels[10]) {
    m_labels = labels;
    refreshDisplay();
}

void ClassificationDisplay::setCurrent(int idx) {
    m_current = idx;
    refreshDisplay();
}

int ClassificationDisplay::current() const {
    return m_current;
}

void ClassificationDisplay::setTypeEdit(QLineEdit *typeEdit) {
    m_typeEdit = typeEdit;
}

void ClassificationDisplay::setSaveCallback(std::function<void()> cb) {
    m_saveCallback = std::move(cb);
}

bool ClassificationDisplay::eventFilter(QObject *watched, QEvent *event) {
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
                if (m_saveCallback) m_saveCallback();
                return true;
            }
        }
    }
    return QObject::eventFilter(watched, event);
}

// ---------------------------------------------------------------------------
// PeriodScrollKeyFilter
// ---------------------------------------------------------------------------

PeriodScrollKeyFilter::PeriodScrollKeyFilter(std::function<void(int)> start, std::function<void()> stop,
                                             std::function<void()> zoomIn, std::function<void()> zoomOut,
                                             QObject *parent)
    : QObject(parent), m_start(std::move(start)), m_stop(std::move(stop)),
      m_zoomIn(std::move(zoomIn)), m_zoomOut(std::move(zoomOut)) {}

bool PeriodScrollKeyFilter::eventFilter(QObject *watched, QEvent *event) {
    if (event->type() != QEvent::KeyPress && event->type() != QEvent::KeyRelease)
        return QObject::eventFilter(watched, event);
    QKeyEvent *ke = static_cast<QKeyEvent *>(event);
    int key = ke->key();
    bool isPlus = (key == Qt::Key_Plus || key == Qt::Key_Equal);
    if (!isPlus && key != Qt::Key_Minus)
        return QObject::eventFilter(watched, event);
    QWidget *focus = QApplication::focusWidget();
    bool textHasFocus = focus && (qobject_cast<QLineEdit *>(focus) || qobject_cast<QAbstractSpinBox *>(focus));
    if (textHasFocus)
        return QObject::eventFilter(watched, event);
    if (ke->modifiers() & Qt::ControlModifier) {
        if (event->type() == QEvent::KeyPress && !ke->isAutoRepeat()) {
            if (isPlus) m_zoomIn(); else m_zoomOut();
        }
        return true;
    }
    if (!ke->isAutoRepeat()) {
        if (event->type() == QEvent::KeyPress)
            m_start(isPlus ? +1 : -1);
        else
            m_stop();
    }
    return true;
}

// ---------------------------------------------------------------------------
// PeriodModifyKeyFilter
// ---------------------------------------------------------------------------

PeriodModifyKeyFilter::PeriodModifyKeyFilter(std::function<void(double)> setFreq,
                                             std::function<void()> toggleDisplay,
                                             std::function<double()> getFreq,
                                             QObject *parent)
    : QObject(parent), m_setFreq(std::move(setFreq)), m_toggleDisplay(std::move(toggleDisplay)), m_getFreq(std::move(getFreq)) {}

bool PeriodModifyKeyFilter::eventFilter(QObject *watched, QEvent *event) {
    if (event->type() != QEvent::KeyPress && event->type() != QEvent::KeyRelease)
        return QObject::eventFilter(watched, event);

    QKeyEvent *ke = static_cast<QKeyEvent *>(event);
    int key = ke->key();

    QWidget *focus = QApplication::focusWidget();
    bool textHasFocus = focus && (qobject_cast<QLineEdit *>(focus) || qobject_cast<QAbstractSpinBox *>(focus));
    if (textHasFocus) return QObject::eventFilter(watched, event);

    if (ke->isAutoRepeat()) {
        if (key == Qt::Key_Asterisk || key == Qt::Key_Slash || key == Qt::Key_Backslash ||
            ((m_multKeyHeld || m_divKeyHeld) && key >= Qt::Key_0 && key <= Qt::Key_9)) {
            return true;
        }
        return QObject::eventFilter(watched, event);
    }

    if (event->type() == QEvent::KeyPress) {
        if (key == Qt::Key_Backslash) {
            m_toggleDisplay();
            return true;
        }
        if (key == Qt::Key_Asterisk && !m_divKeyHeld) {
            m_multKeyHeld = true;
            m_numberPressed = false;
            m_baseFreq = m_getFreq();
            return true;
        }
        if (key == Qt::Key_Slash && !m_multKeyHeld) {
            m_divKeyHeld = true;
            m_numberPressed = false;
            m_baseFreq = m_getFreq();
            return true;
        }
        if ((m_multKeyHeld || m_divKeyHeld) && key >= Qt::Key_1 && key <= Qt::Key_9) {
            m_numberPressed = true;
            double factor = key - Qt::Key_0;
            if (m_multKeyHeld) {
                m_setFreq(m_baseFreq / factor);
            } else if (m_divKeyHeld) {
                m_setFreq(m_baseFreq * factor);
            }
            return true;
        }
        if ((m_multKeyHeld || m_divKeyHeld) && key == Qt::Key_0) {
            return true;
        }
    } else if (event->type() == QEvent::KeyRelease) {
        if (key == Qt::Key_Asterisk && m_multKeyHeld) {
            m_multKeyHeld = false;
            if (!m_numberPressed) m_setFreq(m_baseFreq / 2.0);
            return true;
        }
        if (key == Qt::Key_Slash && m_divKeyHeld) {
            m_divKeyHeld = false;
            if (!m_numberPressed) m_setFreq(m_baseFreq * 2.0);
            return true;
        }
    }

    if ((m_multKeyHeld || m_divKeyHeld) && key >= Qt::Key_0 && key <= Qt::Key_9) {
        return true;
    }
    return QObject::eventFilter(watched, event);
}

// ---------------------------------------------------------------------------
// MainWindow
// ---------------------------------------------------------------------------

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
      m_numpadNav(false),
      m_specFmin(0.0),
      m_specFstep(1.0),
      m_activeMethodIdx(-1),
      m_workerThread(nullptr),
      m_task(nullptr),
      m_modelBusy(false),
      m_modelRefreshPending(false),
      m_pendingModelFreq(-1.0),
      m_stepDir(0) {

    QString defaultLabels[10] = {"nonvar", "var", "unknown", "unknown", "unknown",
                                "unknown", "unknown", "unknown", "unknown", "unknown"};
    for (int i = 0; i < 10; ++i) {
        m_labels[i] = defaultLabels[i];
        m_counts[i] = 0;
    }
    memset(&m_lcData, 0, sizeof(m_lcData));

    loadConfig();

    m_progress = lc_progress_create();
    m_computeCtx = lc_compute_ctx_create(0 /* auto: physical cores */);

    setupUi();
    setupWiring();
}

MainWindow::~MainWindow() {
    if (m_stepTimer) {
        m_stepTimer->stop();
    }
    if (m_progressTimer) {
        m_progressTimer->stop();
    }
    if (m_modelThread) {
        m_modelThread->quit();
        m_modelThread->wait();
    }
    if (m_workerThread) {
        m_workerThread->quit();
        m_workerThread->wait();
    }
    lc_compute_ctx_destroy(m_computeCtx);
    lc_progress_destroy(m_progress);
    lc_free(&m_lcData);
}

void MainWindow::loadConfig() {
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
                    if (s) m_labels[i] = QString::fromUtf8(s->string, s->string_size);
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
                if (!entry.path.isEmpty()) m_files.append(entry);
            }
            continue;
        }

        /* "targets" as JSON array of {path, class} objects */
        if (strcmp(key, "targets") == 0 && el->value->type == json_type_array) {
            json_array_s *arr = json_value_as_array(el->value);
            for (json_array_element_s *ae = arr->start; ae; ae = ae->next) {
                if (ae->value->type != json_type_object) continue;
                json_object_s *tobj = json_value_as_object(ae->value);
                TargetEntry entry;
                for (json_object_element_s *te = tobj->start; te; te = te->next) {
                    if (strcmp(te->name->string, "path") == 0 && te->value->type == json_type_string) {
                        json_string_s *s = json_value_as_string(te->value);
                        if (s) entry.path = QString::fromUtf8(s->string, s->string_size);
                    } else if (strcmp(te->name->string, "class") == 0 && te->value->type == json_type_number) {
                        json_number_s *num = json_value_as_number(te->value);
                        if (num) entry.cls = atoi(num->number);
                    }
                }
                if (!entry.path.isEmpty()) m_targets.append(entry);
            }
            continue;
        }

        /* "period_search" as JSON object */
        if (strcmp(key, "period_search") == 0 && el->value->type == json_type_object) {
            json_object_s *pobj = json_value_as_object(el->value);
            for (json_object_element_s *pe = pobj->start; pe; pe = pe->next) {
                const char *pk = pe->name->string;
                if (pe->value->type == json_type_string && strcmp(pk, "auto_center") == 0) {
                    json_string_s *s = json_value_as_string(pe->value);
                    if (s) m_ps.autoCenter = (strcmp(s->string, "true") == 0);
                    continue;
                }
                if (pe->value->type == json_type_string && strcmp(pk, "display_frequency") == 0) {
                    json_string_s *s = json_value_as_string(pe->value);
                    if (s) m_ps.displayFrequency = (strcmp(s->string, "true") == 0);
                    continue;
                }
                if (pe->value->type == json_type_string && strcmp(pk, "statistic") == 0) {
                    json_string_s *s = json_value_as_string(pe->value);
                    if (s) m_ps.statistic = QString::fromUtf8(s->string, s->string_size);
                    continue;
                }
                if (pe->value->type != json_type_number) continue;
                json_number_s *num = json_value_as_number(pe->value);
                if (!num) continue;
                if (strcmp(pk, "nterms") == 0) m_ps.nterms = atoi(num->number);
                else if (strcmp(pk, "oversampling") == 0) m_ps.oversampling = atof(num->number);
                else if (strcmp(pk, "fmin") == 0) m_ps.fmin = atof(num->number);
                else if (strcmp(pk, "fmax") == 0) m_ps.fmax = atof(num->number);
                else if (strcmp(pk, "zoom_factor") == 0) m_ps.zoomFactor = atof(num->number);
                else if (strcmp(pk, "search_radius") == 0) m_ps.searchRadius = atof(num->number);
                else if (strcmp(pk, "pswf") == 0) m_ps.pswf = atoi(num->number);
                else if (strcmp(pk, "oversmoothing") == 0) m_ps.oversmoothing = atof(num->number);
                else if (strcmp(pk, "nbins") == 0) m_ps.nbins = atoi(num->number);
                else if (strcmp(pk, "scroll_rate") == 0) m_ps.scrollRate = atof(num->number);
            }
            continue;
        }

        if (el->value->type == json_type_string) {
            json_string_s *s = json_value_as_string(el->value);
            if (!s) continue;
            QString val = QString::fromUtf8(s->string, s->string_size);

            if (strcmp(key, "numpad_nav") == 0) {
                m_numpadNav = (val == "true");
            }

            /* Legacy format: "label0"..."label9" */
            for (int i = 0; i < 10; ++i) {
                char lkey[8];
                snprintf(lkey, sizeof(lkey), "label%d", i);
                if (strcmp(key, lkey) == 0) {
                    m_labels[i] = val;
                    break;
                }
            }
        }
    }
    free(root);
}

void MainWindow::saveConfig() {
    QString json = "{\"labels\":[";
    for (int i = 0; i < 10; ++i) {
        json += QString("\"%1\"").arg(escapeJsonString(m_labels[i]));
        if (i < 9) json += ',';
    }
    json += QString("],\"numpad_nav\":\"%1\"").arg(m_numpadNav ? "true" : "false");

    json += ",\"files\":[";
    for (int i = 0; i < m_files.size(); ++i) {
        json += QString("{\"path\":\"%1\",\"n\":%2}").arg(escapeJsonString(m_files[i].path)).arg(m_files[i].n);
        if (i < m_files.size() - 1) json += ',';
    }
    json += "]";

    json += ",\"targets\":[";
    for (int i = 0; i < m_targets.size(); ++i) {
        json += QString("{\"path\":\"%1\",\"class\":%2}").arg(escapeJsonString(m_targets[i].path)).arg(m_targets[i].cls);
        if (i < m_targets.size() - 1) json += ',';
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
                    "\"nbins\":%9,"
                    "\"scroll_rate\":%10,"
                    "\"auto_center\":\"%11\","
                    "\"display_frequency\":\"%12\","
                    "\"statistic\":\"%13\"}")
                .arg(m_ps.nterms)
                .arg(m_ps.oversampling, 0, 'g', 10)
                .arg(m_ps.fmin, 0, 'g', 10)
                .arg(m_ps.fmax, 0, 'g', 10)
                .arg(m_ps.zoomFactor, 0, 'g', 10)
                .arg(m_ps.searchRadius, 0, 'g', 10)
                .arg(m_ps.pswf)
                .arg(m_ps.oversmoothing, 0, 'g', 10)
                .arg(m_ps.nbins)
                .arg(m_ps.scrollRate, 0, 'g', 10)
                .arg(m_ps.autoCenter ? "true" : "false")
                .arg(m_ps.displayFrequency ? "true" : "false")
                .arg(escapeJsonString(m_ps.statistic));
    json += "}";

    QFile f(CONFIG_FILE);
    if (!f.open(QIODevice::WriteOnly))
        return;
    f.write(json.toUtf8());
    f.close();
}

void MainWindow::setupUi() {
    setWindowTitle("lc-qt v ?.?? 25/08/2026");
    resize(900, 850);

    // Menu Bar
    QMenuBar *menuBarPtr = menuBar();
    QMenu *fileMenu = menuBarPtr->addMenu("File");
    m_loadLcAction = fileMenu->addAction("Load Light Curve...");
    fileMenu->addAction("Save Light Curve");
    fileMenu->addAction("Save Light Curve As...");
    fileMenu->addAction("Open Object List...");
    fileMenu->addAction("Save Object List...");
    fileMenu->addSeparator();
    fileMenu->addAction("Quit", qApp, &QCoreApplication::quit);

    QMenu *optionsMenu = menuBarPtr->addMenu("Options");
    optionsMenu->addAction("Customize Period Search...");
    m_periodSearchAction = optionsMenu->actions().last();
    optionsMenu->addAction("Plot Options");
    optionsMenu->addAction("Period Scroll");
    m_periodScrollAction = optionsMenu->actions().last();
    optionsMenu->addAction("Auto Power Spec Calc");

    menuBarPtr->addAction("Help");

    // Main layout
    QWidget *centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);
    mainLayout->setContentsMargins(6, 6, 6, 6);
    mainLayout->setSpacing(6);

    // Global graphite dark style sheet
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
        "QPushButton:checked { background-color: #1d4ed8; border-color: #60a5fa; color: #ffffff; }"
        "QGroupBox { border: 1px solid #32323b; border-radius: 6px; margin-top: 10px; padding: 6px; font-weight: bold; color: #60a5fa; }"
        "QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; background-color: #242429; }"
        "QLabel { color: #cbd5e1; }"
    );

    // Top Bar: Type / HDJ0 / Med,MAD,Amp / Object
    QHBoxLayout *topBarLayout = new QHBoxLayout();
    topBarLayout->addWidget(new QLabel("Type"));
    m_typeEdit = new QLineEdit("unknown");
    m_typeEdit->setMaximumWidth(120);
    m_typeEdit->setReadOnly(true);
    m_typeEdit->setFocusPolicy(Qt::NoFocus);
    topBarLayout->addWidget(m_typeEdit);

    topBarLayout->addWidget(new QLabel("HDJ0"));
    m_hdj0Edit = new QLineEdit("0.000");
    m_hdj0Edit->setMaximumWidth(100);
    m_hdj0Edit->setReadOnly(true);
    m_hdj0Edit->setFocusPolicy(Qt::NoFocus);
    topBarLayout->addWidget(m_hdj0Edit);

    topBarLayout->addWidget(new QLabel("Med"));
    m_medEdit = new QLineEdit("0.000");
    m_medEdit->setMaximumWidth(58);
    m_medEdit->setReadOnly(true);
    m_medEdit->setFocusPolicy(Qt::NoFocus);
    topBarLayout->addWidget(m_medEdit);

    QLabel *madLabel = new QLabel("MAD");
    madLabel->setContentsMargins(0, 0, 0, 0);
    topBarLayout->addWidget(madLabel);
    m_madEdit = new QLineEdit("0.000");
    m_madEdit->setMaximumWidth(48);
    m_madEdit->setReadOnly(true);
    m_madEdit->setFocusPolicy(Qt::NoFocus);
    topBarLayout->addWidget(m_madEdit);

    QLabel *ampLabel = new QLabel("Amp");
    ampLabel->setContentsMargins(0, 0, 0, 0);
    topBarLayout->addWidget(ampLabel);
    m_ampEdit = new QLineEdit("0.000");
    m_ampEdit->setMaximumWidth(48);
    m_ampEdit->setReadOnly(true);
    m_ampEdit->setFocusPolicy(Qt::NoFocus);
    topBarLayout->addWidget(m_ampEdit);

    topBarLayout->addWidget(new QLabel("Object"));
    m_objectEdit = new QLineEdit();
    topBarLayout->addWidget(m_objectEdit);

    mainLayout->addLayout(topBarLayout);

    // Raw and Phased Plots Row
    QHBoxLayout *plotsLayout = new QHBoxLayout();
    m_rawPlot = new LightCurvePlotWidget("Raw Light Curve");
    m_phasedPlot = new PhasedLightCurveWidget("Phased Light Curve");
    plotsLayout->addWidget(m_rawPlot);
    plotsLayout->addWidget(m_phasedPlot);
    mainLayout->addLayout(plotsLayout, 3);

    // Detrending / Light Curve row
    QHBoxLayout *listAndCurveLayout = new QHBoxLayout();
    QGroupBox *statsGroupBox = new QGroupBox("Detrending and prewhitening");
    QHBoxLayout *statsLayout = new QHBoxLayout(statsGroupBox);
    statsLayout->addStretch();
    listAndCurveLayout->addWidget(statsGroupBox, 1);

    QGroupBox *lcGroupBox = new QGroupBox("Light Curve");
    QHBoxLayout *lcLayout = new QHBoxLayout(lcGroupBox);
    m_noPointsLabel = new QLabel("No. points    0");
    m_displayModelBtn = new QPushButton("Display model");
    m_displayModelBtn->setCheckable(true);
    m_displayModelBtn->setChecked(true);
    lcLayout->addWidget(m_noPointsLabel);
    lcLayout->addStretch();
    lcLayout->addWidget(m_displayModelBtn);
    listAndCurveLayout->addWidget(lcGroupBox, 1);

    mainLayout->addLayout(listAndCurveLayout);

    // Log Files Section: Object List + Classification
    QHBoxLayout *logFilesLayout = new QHBoxLayout();

    QGroupBox *objGroupBox = new QGroupBox("Object List");
    QHBoxLayout *objLayout = new QHBoxLayout(objGroupBox);
    m_objListEdit = new QLineEdit();
    m_openObjBtn = new QPushButton("Open");
    QLabel *entryNoLabel = new QLabel("Entry No.");
    m_entryNoEdit = new QLineEdit("0");
    m_entryNoEdit->setMaximumWidth(50);
    objLayout->addWidget(m_objListEdit);
    objLayout->addWidget(m_openObjBtn);
    objLayout->addWidget(entryNoLabel);
    objLayout->addWidget(m_entryNoEdit);
    logFilesLayout->addWidget(objGroupBox, 1);

    QGroupBox *classGroupBox = new QGroupBox("Classification");
    QHBoxLayout *classLayout = new QHBoxLayout(classGroupBox);
    m_customizeBtn = new QPushButton("Customize labels");
    m_statsBtn = new QPushButton("Class stats");
    m_classDisplay = new ClassificationDisplay(m_labels, this);
    m_classDisplay->setTypeEdit(m_typeEdit);

    classLayout->addStretch();
    classLayout->addWidget(m_customizeBtn);
    classLayout->addStretch();
    classLayout->addWidget(m_statsBtn);
    classLayout->addStretch();
    classLayout->addWidget(new QLabel("Cur:"));
    classLayout->addWidget(m_classDisplay);
    classLayout->addStretch();
    logFilesLayout->addWidget(classGroupBox, 1);

    mainLayout->addLayout(logFilesLayout);

    // Period Search / Period Modify section
    QHBoxLayout *periodLayout = new QHBoxLayout();

    QGroupBox *searchGroupBox = new QGroupBox("Period Search");
    QVBoxLayout *searchLayout = new QVBoxLayout(searchGroupBox);
    QHBoxLayout *searchBtns = new QHBoxLayout();
    m_aovBtn = new QPushButton("AoV");
    m_ihsBtn = new QPushButton("IHS");
    m_gbBtn = new QPushButton("GB");
    m_blsBtn = new QPushButton("BLS");
    m_stopBtn = new QPushButton("Stop");
    m_moreBtn = new QPushButton("…");
    m_aovBtn->setCheckable(true);
    m_ihsBtn->setCheckable(true);
    m_gbBtn->setCheckable(true);
    m_blsBtn->setCheckable(true);
    searchBtns->addWidget(m_aovBtn);
    searchBtns->addWidget(m_ihsBtn);
    searchBtns->addWidget(m_gbBtn);
    searchBtns->addWidget(m_blsBtn);
    searchBtns->addWidget(m_stopBtn);
    searchBtns->addWidget(m_moreBtn);
    searchLayout->addLayout(searchBtns);
    m_searchPlot = new SpectrumPlotWidget("Negative Log-Likelihood");
    searchLayout->addWidget(m_searchPlot, 1);
    periodLayout->addWidget(searchGroupBox, 1);

    QGroupBox *modifyGroupBox = new QGroupBox("Period Modify");
    QVBoxLayout *modifyLayout = new QVBoxLayout(modifyGroupBox);
    QHBoxLayout *modifyBtns = new QHBoxLayout();
    m_stepLeftBtn = new QPushButton("<");
    m_periodVal = new QLineEdit("0.000");
    m_periodVal->setMaximumWidth(80);
    m_stepRightBtn = new QPushButton(">");
    m_x2Btn = new QPushButton("x2");
    m_d2Btn = new QPushButton("/2");
    m_otherBtn = new QPushButton("Other");
    m_otherMenu = new QMenu(m_otherBtn);
    m_otherBtn->setMenu(m_otherMenu);
    modifyBtns->addWidget(m_stepLeftBtn);
    modifyBtns->addWidget(m_periodVal);
    modifyBtns->addWidget(m_stepRightBtn);
    modifyBtns->addWidget(m_x2Btn);
    modifyBtns->addWidget(m_d2Btn);
    modifyBtns->addWidget(m_otherBtn);
    modifyBtns->addWidget(new QPushButton("◄"));
    modifyBtns->addWidget(new QPushButton("►"));
    modifyLayout->addLayout(modifyBtns);

    QHBoxLayout *bottomModifyArea = new QHBoxLayout();
    m_minusBtn = new QPushButton("-");
    m_minusBtn->setMaximumSize(25, 50);
    m_plusBtn = new QPushButton("+");
    m_plusBtn->setMaximumSize(25, 50);
    m_zoomPlot = new ZoomedSpectrumWidget("Zoomed Spectrum");

    QVBoxLayout *leftCtrl = new QVBoxLayout();
    leftCtrl->addWidget(m_minusBtn);
    leftCtrl->addStretch();

    QVBoxLayout *rightCtrl = new QVBoxLayout();
    rightCtrl->addWidget(m_plusBtn);
    rightCtrl->addStretch();

    bottomModifyArea->addLayout(leftCtrl);
    bottomModifyArea->addWidget(m_zoomPlot, 1);
    bottomModifyArea->addLayout(rightCtrl);

    modifyLayout->addLayout(bottomModifyArea, 1);
    periodLayout->addWidget(modifyGroupBox, 1);

    mainLayout->addLayout(periodLayout, 2);

    // Status bar widgets
    m_phasedTimeLabel = new QLabel();
    m_phasedTimeLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    statusBar()->addPermanentWidget(m_phasedTimeLabel);

    m_progressLabel = new QLabel();
    m_progressLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    statusBar()->addPermanentWidget(m_progressLabel);

    setButtonsEnabled(true);
}

void MainWindow::setupWiring() {
    qApp->installEventFilter(m_classDisplay);

    m_classDisplay->setSaveCallback([this]() {
        if (m_currentPath.isEmpty()) return;
        int cls = m_classDisplay->current();
        bool found = false;
        for (auto &t : m_targets) {
            if (t.path == m_currentPath) { t.cls = cls; found = true; break; }
        }
        if (!found) m_targets.append(TargetEntry{m_currentPath, cls});
        saveConfig();
    });

    connect(m_loadLcAction, &QAction::triggered, this, [this]() {
        QString path = QFileDialog::getOpenFileName(
            this, "Load Light Curve", QString(),
            "Light Curve Files (*.dat *.fits);;All Files (*)");
        loadLightCurve(path);
    });

    connect(m_customizeBtn, &QPushButton::clicked, this, [this]() {
        bool reopen = true;
        while (reopen) {
            reopen = false;
            CustomizeLabelsDialog dlg(m_labels, &m_numpadNav, this);
            connect(&dlg, &CustomizeLabelsDialog::labelsChanged, m_classDisplay, [this]() {
                m_classDisplay->refreshDisplay();
                saveConfig();
            });
            dlg.exec();
            if (dlg.shouldReopen()) {
                reopen = true;
            }
        }
    });

    connect(m_statsBtn, &QPushButton::clicked, this, [this]() {
        ClassificationStatsDialog dlg(m_labels, m_counts, this);
        dlg.exec();
    });

    // Zoomed spectrum wiring
    connect(m_minusBtn, &QPushButton::clicked, m_zoomPlot, &ZoomedSpectrumWidget::zoomOut);
    connect(m_plusBtn, &QPushButton::clicked, m_zoomPlot, &ZoomedSpectrumWidget::zoomIn);
    connect(m_searchPlot, &SpectrumPlotWidget::frequencyClicked, m_zoomPlot,
            [this](double freq) { m_zoomPlot->selectFrequency(freq, true); });
    connect(m_zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged,
            m_searchPlot, &SpectrumPlotWidget::setSelectedFrequency);
    connect(m_zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged,
            m_phasedPlot, &PhasedLightCurveWidget::setFrequency);
    connect(m_searchPlot, &SpectrumPlotWidget::rangeSelected,
            m_zoomPlot, &ZoomedSpectrumWidget::setViewFromSelection);

    static const struct {
        const char *label;
        double factor;
    } otherMults[] = {{"x3", 1.0 / 3.0}, {"/3", 3.0}, {"x5", 1.0 / 5.0}, {"/5", 5.0}, {"x7", 1.0 / 7.0}, {"/7", 7.0}};
    for (const auto &m : otherMults) {
        QAction *act = m_otherMenu->addAction(QString::fromLatin1(m.label));
        connect(act, &QAction::triggered, this, [this, m] { m_zoomPlot->multiplyFrequency(m.factor); });
    }
    connect(m_x2Btn, &QPushButton::clicked, this, [this] { m_zoomPlot->multiplyFrequency(0.5); });
    connect(m_d2Btn, &QPushButton::clicked, this, [this] { m_zoomPlot->multiplyFrequency(2.0); });

    connect(m_zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged, m_periodVal,
            [this](double freq) {
                double dT = (m_lcData.n > 1) ? m_lcData.x[m_lcData.n - 1] - m_lcData.x[0] : 0.0;
                m_periodVal->setText(formatCoordinateUpstream(freq, dT, m_ps.oversampling, m_ps.nterms, m_ps.displayFrequency));
            });
    connect(m_periodVal, &QLineEdit::editingFinished, this, [this]() {
        bool ok = false;
        double val = m_periodVal->text().toDouble(&ok);
        if (ok && val > 0.0) {
            double f = m_ps.displayFrequency ? val : (1.0 / val);
            m_zoomPlot->selectFrequency(f, false);
        }
    });

    for (QLineEdit *field : {m_objectEdit, m_periodVal, m_entryNoEdit, m_objListEdit}) {
        connect(field, &QLineEdit::editingFinished, this, [field]() {
            field->deselect();
            field->clearFocus();
        });
    }

    m_stepTimer = new QTimer(this);
    m_stepTimer->setInterval(20);
    connect(m_stepTimer, &QTimer::timeout, this, [this]() {
        double dt = (double)m_stepClock.restart() * 1e-3;
        double dT = (m_lcData.n > 1) ? m_lcData.x[m_lcData.n - 1] - m_lcData.x[0] : 0.0;
        if (dT > 0.0 && m_stepDir != 0) m_zoomPlot->stepFrequency((double)m_stepDir * m_ps.scrollRate * dt / dT);
    });

    connect(m_stepLeftBtn, &QPushButton::pressed, this, [this] { startStep(-1); });
    connect(m_stepLeftBtn, &QPushButton::released, this, [this] { stopStep(); });
    connect(m_stepRightBtn, &QPushButton::pressed, this, [this] { startStep(+1); });
    connect(m_stepRightBtn, &QPushButton::released, this, [this] { stopStep(); });

    m_scrollFilter = new PeriodScrollKeyFilter([this](int dir) { startStep(dir); },
                                               [this]() { stopStep(); },
                                               [this]() { m_zoomPlot->zoomIn(); },
                                               [this]() { m_zoomPlot->zoomOut(); },
                                               this);
    qApp->installEventFilter(m_scrollFilter);

    m_modifyFilter = new PeriodModifyKeyFilter(
        [this](double f) { m_zoomPlot->selectFrequency(f, false); },
        [this]() {
            m_ps.displayFrequency = !m_ps.displayFrequency;
            saveConfig();
            double freq = m_zoomPlot->centerFrequency();
            double dT = (m_lcData.n > 1) ? m_lcData.x[m_lcData.n - 1] - m_lcData.x[0] : 0.0;
            if (freq > 0.0) {
                m_periodVal->setText(formatCoordinateUpstream(freq, dT, m_ps.oversampling, m_ps.nterms, m_ps.displayFrequency));
            }
        },
        [this]() { return m_zoomPlot->centerFrequency(); },
        this
    );
    qApp->installEventFilter(m_modifyFilter);

    // Phased model worker
    m_modelThread = new QThread(this);
    m_modelWorker = new PhasedModelWorker();
    m_modelWorker->moveToThread(m_modelThread);
    connect(m_modelThread, &QThread::finished, m_modelWorker, &QObject::deleteLater);
    m_modelThread->start();

    connect(m_modelWorker, &PhasedModelWorker::done, this,
            [this](int rc, double freq, QVector<float> model, int style, qint64 elapsedNs) {
                m_modelBusy = false;
                m_phasedTimeLabel->setText(QString("Phased model: %1 ms").arg((double)elapsedNs * 1e-6, 0, 'f', 2));
                if (rc == 0 && m_lcData.n > 0 && model.size() == (int)m_lcData.n) {
                    m_phasedPlot->setModel(model.constData(), m_lcData.n, (lc_model_style_t)style, freq);
                    float mn = model[0], mx = model[0];
                    for (int i = 1; i < model.size(); ++i) {
                        if (model[i] < mn) mn = model[i];
                        if (model[i] > mx) mx = model[i];
                    }
                    m_ampEdit->setText(QString::number((double)(mx - mn), 'f', 3));
                }
                tryDispatchModel();
            });

    connect(m_zoomPlot, &ZoomedSpectrumWidget::centerFrequencyChanged, this, [this](double freq) {
        if (m_lcData.n == 0 || m_activeMethodIdx < 0 || !(freq > 0.0)) {
            m_phasedPlot->clearModel();
            m_ampEdit->setText(QStringLiteral("0.000"));
            m_modelRefreshPending = false;
            return;
        }
        m_pendingModelFreq = freq;
        m_modelRefreshPending = true;
        tryDispatchModel();
    });

    connect(m_displayModelBtn, &QPushButton::toggled, m_phasedPlot, &PhasedLightCurveWidget::setDisplayModel);

    // Progress timer
    m_progressTimer = new QTimer(this);
    m_progressTimer->setInterval(100);
    connect(m_progressTimer, &QTimer::timeout, this, [this]() {
        uint32_t done = lc_progress_done(m_progress);
        uint32_t total = lc_progress_total(m_progress);
        if (total == 0) {
            m_searchPlot->setProgressText("...");
            m_progressLabel->setText("Computation in progress...");
            return;
        }
        double pct = 100.0 * (double)done / (double)total;
        m_searchPlot->setProgressText(QString("%1%").arg(pct, 0, 'f', 1));
        QString timeLeft;
        double elapsed = (double)m_elapsedTimer.elapsed() / 1000.0;
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
        m_progressLabel->setText(QString("Computation in progress: %1% complete | %2").arg(pct, 0, 'f', 1).arg(timeLeft));
    });

    connect(m_aovBtn, &QPushButton::clicked, this, [this]() { launchSpectrum(LC_SPEC_AOV); });
    connect(m_ihsBtn, &QPushButton::clicked, this, [this]() { launchSpectrum(LC_SPEC_IHS); });
    connect(m_gbBtn, &QPushButton::clicked, this, [this]() { launchSpectrum(LC_SPEC_GB); });
    connect(m_blsBtn, &QPushButton::clicked, this, [this]() { launchSpectrum(LC_SPEC_BLS); });
    connect(m_stopBtn, &QPushButton::clicked, this, [this]() { lc_progress_request_cancel(m_progress); });

    // MMB peak search wiring
    connect(m_searchPlot, &SpectrumPlotWidget::middleClicked, this, [this](double freq) {
        if (m_specNll.isEmpty()) return;
        double fullSpan = (m_specNll.size() > 1) ? (double)(m_specNll.size() - 1) * m_specFstep : 1.0;
        double halfWidth = m_ps.searchRadius * fullSpan;
        double peakFreq = findPeakNear(freq, halfWidth);
        selectPeakAndAlignExtremum(peakFreq);
    });

    connect(m_zoomPlot, &ZoomedSpectrumWidget::middleClicked, this, [this](double freq, double fov) {
        if (m_specNll.isEmpty()) return;
        double halfWidth = fov / 2.0;
        double peakFreq = findPeakNear(freq, halfWidth);
        selectPeakAndAlignExtremum(peakFreq);
    });

    connect(m_moreBtn, &QPushButton::clicked, this, &MainWindow::openPeriodSearchDialog);
    connect(m_periodSearchAction, &QAction::triggered, this, &MainWindow::openPeriodSearchDialog);
    connect(m_periodScrollAction, &QAction::triggered, this, &MainWindow::openPeriodScrollDialog);
}

void MainWindow::startStep(int dir) {
    m_stepDir = dir;
    m_stepClock.start();
    m_stepTimer->start();
}

void MainWindow::stopStep() {
    m_stepTimer->stop();
    m_stepDir = 0;
}

void MainWindow::setButtonsEnabled(bool enabled) {
    bool hasData = (m_lcData.n > 0);
    m_aovBtn->setEnabled(enabled && hasData);
    m_ihsBtn->setEnabled(enabled && hasData);
    m_gbBtn->setEnabled(enabled && hasData);
    m_blsBtn->setEnabled(enabled && hasData);
    m_stopBtn->setEnabled(!enabled);
}

void MainWindow::loadLightCurve(const QString &path) {
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
        QMessageBox::warning(this, "Load Error",
                             QString("Failed to load light curve from:\n%1").arg(path));
        return;
    }

    lc_detrend(&newData);

    lc_free(&m_lcData);
    m_lcData = newData;

    m_rawPlot->setData(&m_lcData);
    m_phasedPlot->setData(&m_lcData);

    m_noPointsLabel->setText(QString("No. points    %1").arg(m_lcData.n));
    m_objectEdit->setText(fi.fileName());
    m_currentPath = path;

    m_hdj0Edit->setText(QString::number(m_lcData.x[0], 'f', 3));

    QVector<float> sorted(m_lcData.y, m_lcData.y + m_lcData.n);
    std::sort(sorted.begin(), sorted.end());
    double median;
    if (sorted.size() % 2 == 0)
        median = ((double)sorted[sorted.size() / 2 - 1] + (double)sorted[sorted.size() / 2]) / 2.0;
    else
        median = (double)sorted[sorted.size() / 2];
    m_medEdit->setText(QString::number(median, 'f', 3));

    double mad = 0.0;
    for (unsigned int i = 0; i < m_lcData.n; i++)
        mad += std::fabs((double)m_lcData.y[i] - median);
    mad /= (double)m_lcData.n;
    m_madEdit->setText(QString::number(mad, 'f', 3));

    bool found = false;
    for (int i = 0; i < m_files.size(); i++) {
        if (m_files[i].path == path) {
            m_files[i].n = (int)m_lcData.n;
            found = true;
            break;
        }
    }
    if (!found) {
        m_files.append(FileEntry{path, (int)m_lcData.n});
    }
    saveConfig();
    setButtonsEnabled(true);
}

void MainWindow::tryDispatchModel() {
    if (m_modelBusy) return;
    if (!m_modelRefreshPending) return;
    if (m_lcData.n == 0 || m_activeMethodIdx < 0 || !(m_pendingModelFreq > 0.0)) {
        m_modelRefreshPending = false;
        return;
    }
    m_modelRefreshPending = false;
    m_modelBusy = true;

    lc_periodogram_config_t mcfg;
    memset(&mcfg, 0, sizeof(mcfg));
    mcfg.method = (lc_spec_method_t)m_activeMethodIdx;
    mcfg.nterms = m_ps.nterms;
    mcfg.oversampling = m_ps.oversampling;
    mcfg.fmin = m_ps.fmin;
    mcfg.fmax = m_ps.fmax;
    mcfg.pswf = m_ps.pswf;
    mcfg.oversmoothing = m_ps.oversmoothing;
    mcfg.nbins = m_ps.nbins;

    lc_data_t clone = cloneLcData(&m_lcData);
    double freq = m_pendingModelFreq;
    QMetaObject::invokeMethod(m_modelWorker, [this, clone, mcfg, freq]() {
        m_modelWorker->compute(clone, mcfg, freq);
    }, Qt::QueuedConnection);
}

void MainWindow::moveToPeak(double peakFreq) {
    double dT = (m_lcData.n > 1) ? m_lcData.x[m_lcData.n - 1] - m_lcData.x[0] : 0.0;
    if (dT > 0.0) {
        double halfFov = 2.0 / dT;
        m_zoomPlot->setViewFromSelection(peakFreq - halfFov, peakFreq + halfFov);
    } else {
        m_zoomPlot->selectFrequency(peakFreq, false);
    }
}

void MainWindow::selectPeakAndAlignExtremum(double peakFreq) {
    if (m_activeMethodIdx >= 0 && m_lcData.n > 0 && peakFreq > 0.0) {
        lc_periodogram_config_t mcfg;
        memset(&mcfg, 0, sizeof(mcfg));
        mcfg.method = (lc_spec_method_t)m_activeMethodIdx;
        mcfg.nterms = m_ps.nterms;
        mcfg.oversampling = m_ps.oversampling;
        mcfg.fmin = m_ps.fmin;
        mcfg.fmax = m_ps.fmax;
        mcfg.pswf = m_ps.pswf;
        mcfg.oversmoothing = m_ps.oversmoothing;
        mcfg.nbins = m_ps.nbins;

        QVector<float> model(m_lcData.n);
        lc_model_style_t style = LC_MODEL_LINE;
        int rc = lc_compute_phased_model(&m_lcData, &mcfg, peakFreq, model.data(), &style);
        if (rc == 0) {
            double offset = lc_compute_phase_offset(&m_lcData, &mcfg, peakFreq, model.constData());
            m_phasedPlot->setPhaseOffset(offset);
            m_phasedPlot->setModel(model.constData(), m_lcData.n, style, peakFreq);
            float mn = model[0], mx = model[0];
            for (int i = 1; i < model.size(); ++i) {
                if (model[i] < mn) mn = model[i];
                if (model[i] > mx) mx = model[i];
            }
            m_ampEdit->setText(QString::number((double)(mx - mn), 'f', 3));
        }
    }
    moveToPeak(peakFreq);
}

void MainWindow::launchSpectrum(lc_spec_method_t method) {
    if (m_lcData.n == 0) {
        QMessageBox::warning(this, "No Data", "Load a light curve first.");
        return;
    }
    if (m_workerThread) return;

    lc_periodogram_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.method = method;
    cfg.statistic = (m_ps.statistic == "raw") ? LC_STAT_RAW : LC_STAT_BAYES;
    cfg.nterms = m_ps.nterms;
    cfg.oversampling = m_ps.oversampling;
    cfg.fmin = m_ps.fmin;
    cfg.fmax = m_ps.fmax;
    cfg.pswf = m_ps.pswf;
    cfg.nthreads = 0;
    cfg.oversmoothing = m_ps.oversmoothing;
    cfg.nbins = m_ps.nbins;
    cfg.peak_threshold = 8.0;

    std::array<QPushButton*, 4> methodBtns = {m_ihsBtn, m_aovBtn, m_gbBtn, m_blsBtn};
    int methodIdx = (int)method;
    int prevMethodIdx = m_activeMethodIdx;
    for (int i = 0; i < 4; ++i) methodBtns[i]->setChecked(i == methodIdx);

    lc_progress_reset(m_progress);
    setButtonsEnabled(false);
    m_elapsedTimer.start();
    m_searchPlot->setProgressText("...");
    m_progressLabel->setText("Computation in progress...");
    m_progressTimer->start();

    m_workerThread = new QThread();
    m_task = new PeriodogramTask(m_computeCtx, &m_lcData, cfg, m_progress);
    m_task->moveToThread(m_workerThread);

    connect(m_workerThread, &QThread::started, m_task, &PeriodogramTask::run);
    connect(m_task, &PeriodogramTask::finished, this,
            [this, methodIdx](double fmin, double fstep, QVector<float> nll, QVector<QPair<double, float>>) {
                m_activeMethodIdx = methodIdx;
                m_searchPlot->setData(fmin, fstep, nll);
                m_specFmin = fmin;
                m_specFstep = fstep;
                m_specNll = nll;
            });
    connect(m_task, &PeriodogramTask::finished, m_zoomPlot,
            [this](double fmin, double fstep, QVector<float> nll, QVector<QPair<double, float>>) {
                m_zoomPlot->setZoomFactor(m_ps.zoomFactor);
                m_zoomPlot->setFullSpectrum(fmin, fstep, nll);
            });
    connect(m_task, &PeriodogramTask::finished, this,
            [this](double, double, QVector<float>, QVector<QPair<double, float>> peaks) {
                m_detectedPeaks = peaks;
                if (m_ps.autoCenter && !peaks.isEmpty()) {
                    selectPeakAndAlignExtremum(peaks[0].first);
                }
            });
    connect(m_task, &PeriodogramTask::failed, this,
            [this](QString msg) { QMessageBox::warning(this, "Computation Failed", msg); });

    auto onComplete = [this]() {
        m_progressTimer->stop();
        m_progressLabel->clear();
        m_searchPlot->setProgressText("");
        setButtonsEnabled(true);
        if (m_workerThread) m_workerThread->quit();
    };
    connect(m_task, &PeriodogramTask::finished, this, [onComplete](double, double, QVector<float>, QVector<QPair<double, float>>) {
        onComplete();
    });
    connect(m_task, &PeriodogramTask::cancelled, this, [this, onComplete, prevMethodIdx, methodBtns]() {
        for (int i = 0; i < 4; ++i) methodBtns[i]->setChecked(i == prevMethodIdx);
        onComplete();
    });
    connect(m_task, &PeriodogramTask::failed, this, [this, onComplete, prevMethodIdx, methodBtns](QString) {
        for (int i = 0; i < 4; ++i) methodBtns[i]->setChecked(i == prevMethodIdx);
        onComplete();
    });

    connect(m_workerThread, &QThread::finished, m_task, &QObject::deleteLater);
    connect(m_workerThread, &QThread::finished, m_workerThread, &QObject::deleteLater);
    connect(m_workerThread, &QThread::finished, this, [this]() {
        m_workerThread = nullptr;
        m_task = nullptr;
    });

    m_workerThread->start();
}

static void quadraticVertex(double freq, double fstep, float left, float center, float right, double oversampling, double *peakFreq, float *peakNll) {
    *peakFreq = freq;
    *peakNll = center;
    if (oversampling < 3.0) return;
    double x0 = freq - fstep, x1 = freq, x2 = freq + fstep;
    double slope01 = ((double)center - (double)left) / (x1 - x0);
    double slope12 = ((double)right - (double)center) / (x2 - x1);
    double curvature = (slope12 - slope01) / (x2 - x0);
    if (curvature == 0.0) return;
    double linear = slope01 - curvature * (x0 + x1);
    double vx = -linear / (2.0 * curvature);
    double vy = (double)left + slope01 * (vx - x0) + curvature * (vx - x0) * (vx - x1);
    if (std::isfinite(vx) && std::isfinite(vy) && vx >= x0 && vx <= x2) {
        *peakFreq = vx;
        *peakNll = (float)vy;
    }
}

double MainWindow::findPeakNear(double clickFreq, double searchHalfWidth) {
    if (m_specNll.isEmpty()) return clickFreq;
    double lo = clickFreq - searchHalfWidth;
    double hi = clickFreq + searchHalfWidth;

    for (const auto &pk : m_detectedPeaks) {
        if (pk.first >= lo && pk.first <= hi) return pk.first;
    }

    int i0 = (int)std::ceil((lo - m_specFmin) / m_specFstep);
    int i1 = (int)std::floor((hi - m_specFmin) / m_specFstep);
    if (i0 < 1) i0 = 1;
    if (i1 > m_specNll.size() - 2) i1 = m_specNll.size() - 2;

    double bestFreq = clickFreq;
    float bestNll = -1.0f;
    bool found = false;
    double oversampling = m_ps.oversampling;

    for (int i = i0; i <= i1; ++i) {
        float left = m_specNll[i - 1], center = m_specNll[i], right = m_specNll[i + 1];
        if (center > left && center > right) {
            double freq = m_specFmin + (double)i * m_specFstep;
            double pf = freq;
            float pn = center;
            quadraticVertex(freq, m_specFstep, left, center, right, oversampling, &pf, &pn);
            if (!found || pn > bestNll) {
                bestFreq = pf;
                bestNll = pn;
                found = true;
            }
        }
    }

    if (!found) {
        for (int i = i0; i <= i1; ++i) {
            if (m_specNll[i] > bestNll) {
                bestNll = m_specNll[i];
                bestFreq = m_specFmin + (double)i * m_specFstep;
                found = true;
            }
        }
    }
    return bestFreq;
}

void MainWindow::openPeriodSearchDialog() {
    CustomizePeriodSearchDialog dlg(&m_ps, this);
    if (dlg.exec() == QDialog::Accepted) {
        saveConfig();
    }
}

void MainWindow::openPeriodScrollDialog() {
    PeriodScrollDialog dlg(&m_ps, this);
    if (dlg.exec() == QDialog::Accepted) {
        saveConfig();
    }
}
