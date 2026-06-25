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
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QDir>

#include "json.h"

#include "windows/customize_labels.h"
#include "windows/classification_stats.h"

static const char *CONFIG_FILE = ".lc-config.json";

static void loadConfig(QString labels[10], bool *numpadNav) {
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
        if (el->value->type == json_type_string) {
            json_string_s *s = json_value_as_string(el->value);
            if (!s) continue;
            QString val = QString::fromUtf8(s->string, s->string_size);

            if (strcmp(el->name->string, "numpad_nav") == 0) {
                *numpadNav = (val == "true");
            }

            for (int i = 0; i < 10; ++i) {
                char key[8];
                snprintf(key, sizeof(key), "label%d", i);
                if (strcmp(el->name->string, key) == 0) {
                    labels[i] = val;
                    break;
                }
            }
        }
    }
    free(root);
}

static void saveConfig(const QString labels[10], bool numpadNav) {
    QString json = "{";
    for (int i = 0; i < 10; ++i) {
        QString escaped = labels[i];
        escaped.replace('\\', "\\\\").replace('"', "\\\"");
        json += QString("\"label%1\":\"%2\"").arg(i).arg(escaped);
        if (i < 9) json += ',';
    }
    json += QString(",\"numpad_nav\":\"%1\"}").arg(numpadNav ? "true" : "false");

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
public:
    ClassificationDisplay(QString labels[10], QWidget *parent = nullptr)
        : QLineEdit(parent), m_labels(labels), m_current(0) {
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

    bool eventFilter(QObject *watched, QEvent *event) override {
        if (event->type() == QEvent::KeyPress) {
            QKeyEvent *ke = static_cast<QKeyEvent *>(event);
            int key = ke->key();
            if (key >= Qt::Key_0 && key <= Qt::Key_9) {
                QWidget *focus = QApplication::focusWidget();
                bool textHasFocus = focus && qobject_cast<QLineEdit *>(focus);
                if (!textHasFocus) {
                    m_current = key - Qt::Key_0;
                    refreshDisplay();
                    return true;
                }
            }
        }
        return QObject::eventFilter(watched, event);
    }
};


int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    
    // Set classic Motif-like color palette and style hint
    app.setStyle("Fusion");
    
    QMainWindow window;
    window.setWindowTitle("lc v ?.?? 25/05/2026");
    window.resize(900, 850);

    // Classification state
    QString labels[10] = {"nonvar", "variable", "variable", "variable", "variable",
                          "variable", "variable", "variable", "variable", "variable"};
    bool numpadNav = false;
    loadConfig(labels, &numpadNav);
    int counts[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    // Menu Bar
    QMenuBar *menuBar = window.menuBar();
    QMenu *fileMenu = menuBar->addMenu("File");
    fileMenu->addAction("Load Light Curve...");
    fileMenu->addAction("Save Light Curve");
    fileMenu->addAction("Save Light Curve As...");
    fileMenu->addAction("Open Object List...");
    fileMenu->addAction("Save Object List...");
    fileMenu->addSeparator();
    fileMenu->addAction("Quit", &app, &QCoreApplication::quit);
    
    QMenu *optionsMenu = menuBar->addMenu("Options");
    optionsMenu->addAction("Customize In/Out...");
    optionsMenu->addAction("Preview");
    optionsMenu->addAction("Classes");
    optionsMenu->addAction("Frequency analysis");
    optionsMenu->addAction("Customize Period Search...");
    optionsMenu->addAction("Light Curve");
    optionsMenu->addAction("Plot Options");
    optionsMenu->addAction("Period Scroll");
    optionsMenu->addAction("Multi-Period Window");
    optionsMenu->addAction("Auto Equalize");
    optionsMenu->addAction("Auto Power Spec Load");
    optionsMenu->addAction("Auto Power Spec Calc");
    optionsMenu->addAction("Auto Adjust Spectrum Max");
    
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
    QLineEdit *typeEdit = new QLineEdit("Var");
    typeEdit->setMaximumWidth(120);
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
    MockPlotWidget *rawPlot = new MockPlotWidget("Raw Light Curve (Scatter & Error Bars)", false, false, false);
    MockPlotWidget *phasedPlot = new MockPlotWidget("Phased Light Curve (Sine Model + Scatter)", true, false, false);
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

    classLayout->addStretch();
    classLayout->addWidget(customizeBtn);
    classLayout->addStretch();
    classLayout->addWidget(new QLabel("Current:"));
    classLayout->addWidget(classDisplay);
    classLayout->addStretch();
    logFilesLayout->addWidget(classGroupBox, 1);

    mainLayout->addLayout(logFilesLayout);

    app.installEventFilter(classDisplay);

    QObject::connect(customizeBtn, &QPushButton::clicked, [&]() {
        CustomizeLabelsDialog dlg(labels, &numpadNav, &window);
        QObject::connect(&dlg, &CustomizeLabelsDialog::labelsChanged, classDisplay, [classDisplay, labels, &numpadNav]() {
            classDisplay->refreshDisplay();
            saveConfig(labels, numpadNav);
        });
        dlg.exec();
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
    searchBtns->addWidget(new QPushButton("AoV"));
    searchBtns->addWidget(new QPushButton("IHS"));
    searchBtns->addWidget(new QPushButton("GB"));
    searchBtns->addWidget(new QPushButton("BLS"));
    searchBtns->addWidget(new QPushButton("P"));
    searchBtns->addWidget(new QPushButton("Stop"));
    searchBtns->addWidget(new QPushButton("Load"));
    searchBtns->addWidget(new QPushButton("..."));
    searchLayout->addLayout(searchBtns);
    MockPlotWidget *searchPlot = new MockPlotWidget("Period Search Spectrum (Power vs Frequency)", false, true, false);
    searchLayout->addWidget(searchPlot, 1);
    periodLayout->addWidget(searchGroupBox, 1);

    
    QGroupBox *modifyGroupBox = new QGroupBox("Period Modify");
    QVBoxLayout *modifyLayout = new QVBoxLayout(modifyGroupBox);
    QHBoxLayout *modifyBtns = new QHBoxLayout();
    QLineEdit *periodVal = new QLineEdit("0.000");
    periodVal->setMaximumWidth(80);
    modifyBtns->addWidget(periodVal);
    modifyBtns->addWidget(new QPushButton("<"));
    modifyBtns->addWidget(new QPushButton(">"));
    modifyBtns->addWidget(new QPushButton("x2"));
    modifyBtns->addWidget(new QPushButton("/2"));
    modifyBtns->addWidget(new QPushButton("x3"));
    modifyBtns->addWidget(new QPushButton("/3"));
    modifyBtns->addWidget(new QPushButton("◄"));
    modifyBtns->addWidget(new QPushButton("►"));
    modifyLayout->addLayout(modifyBtns);
    
    // Bottom layout for modify showing [-] and [+] on the side
    QHBoxLayout *bottomModifyArea = new QHBoxLayout();
    QPushButton *minusBtn = new QPushButton("-");
    minusBtn->setMaximumSize(20, 20);
    QPushButton *plusBtn = new QPushButton("+");
    plusBtn->setMaximumSize(20, 20);
    
    MockPlotWidget *modifyPlot = new MockPlotWidget("Period Modify (Model Fitting)", false, false, true);

    
    QVBoxLayout *leftCtrl = new QVBoxLayout();
    leftCtrl->addWidget(minusBtn);
    leftCtrl->addStretch();
    
    QVBoxLayout *rightCtrl = new QVBoxLayout();
    rightCtrl->addWidget(plusBtn);
    rightCtrl->addStretch();
    
    bottomModifyArea->addLayout(leftCtrl);
    bottomModifyArea->addWidget(modifyPlot, 1);
    bottomModifyArea->addLayout(rightCtrl);
    
    modifyLayout->addLayout(bottomModifyArea, 1);
    
    periodLayout->addWidget(modifyGroupBox, 1);
    
    mainLayout->addLayout(periodLayout, 2); // stretch factor 2
    
    window.show();
    return app.exec();
}
