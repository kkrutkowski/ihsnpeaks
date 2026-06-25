#include "customize_labels.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>

static const int NUMPAD_SKIP[] = {2, 4, 6, 8};

static bool isNumPadSkip(int i) {
    for (int skip : NUMPAD_SKIP) {
        if (i == skip) return true;
    }
    return false;
}

CustomizeLabelsDialog::CustomizeLabelsDialog(QString labels[10], bool *numpadNav, QWidget *parent)
    : QDialog(parent), m_labels(labels), m_numpadNav(numpadNav), m_grid(nullptr), m_shouldReopen(false) {

    setWindowTitle("Customize Labels");
    setMinimumWidth(340);
    setStyleSheet(
        "QDialog { background-color: #242429; color: #e2e8f0; }"
        "QLabel { color: #cbd5e1; }"
        "QLineEdit { background-color: #17171b; border: 1px solid #32323b; border-radius: 4px; "
        "  color: #ffffff; padding: 3px; }"
        "QLineEdit:focus { border: 1px solid #3b82f6; }"
        "QPushButton { background-color: #303038; border: 1px solid #3d3d47; border-radius: 4px; "
        "  color: #f1f5f9; padding: 4px 12px; font-weight: bold; }"
        "QPushButton:hover { background-color: #3b82f6; border-color: #60a5fa; }"
        "QPushButton:pressed { background-color: #1d4ed8; }"
        "QCheckBox { color: #cbd5e1; }"
        "QGroupBox { border: 1px solid #32323b; border-radius: 6px; margin-top: 10px; padding: 6px; "
        "  font-weight: bold; color: #60a5fa; }"
        "QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; "
        "  background-color: #242429; }"
    );

    m_mainLayout = new QVBoxLayout(this);
    m_mainLayout->setContentsMargins(12, 12, 12, 12);
    m_mainLayout->setSpacing(8);

    QHBoxLayout *topLayout = new QHBoxLayout();
    topLayout->addStretch();
    m_numpadCheck = new QCheckBox("NumPad navigation");
    m_numpadCheck->setChecked(*m_numpadNav);
    topLayout->addWidget(m_numpadCheck);
    topLayout->addStretch();
    m_mainLayout->addLayout(topLayout);

    m_gridBox = new QGroupBox("Label Assignments");
    m_mainLayout->addWidget(m_gridBox);

    m_grid = new QGridLayout(m_gridBox);
    m_grid->setSpacing(6);
    m_grid->addWidget(new QLabel("#"), 0, 0);
    m_grid->addWidget(new QLabel("Label"), 0, 1);

    static const char *defaults[] = {
        "nonvar", "var", "unknown", "unknown", "unknown",
        "unknown", "unknown", "unknown", "unknown", "unknown"
    };

    for (int i = 0; i < 10; ++i) {
        m_numLabels[i] = new QLabel(QString::number(i));
        m_numLabels[i]->setMinimumWidth(20);
        m_numLabels[i]->setAlignment(Qt::AlignCenter);

        m_edits[i] = new QLineEdit(m_labels[i].isEmpty() ? QString(defaults[i]) : m_labels[i]);

        m_grid->addWidget(m_numLabels[i], i + 1, 0);
        m_grid->addWidget(m_edits[i], i + 1, 1);
    }

    m_grid->setColumnStretch(1, 1);

    updateVisibility();

    QHBoxLayout *btnLayout = new QHBoxLayout();
    btnLayout->addStretch();
    QPushButton *okBtn = new QPushButton("OK");
    QPushButton *cancelBtn = new QPushButton("Cancel");
    btnLayout->addWidget(okBtn);
    btnLayout->addWidget(cancelBtn);
    m_mainLayout->addLayout(btnLayout);

    connect(m_numpadCheck, &QCheckBox::stateChanged, this, &CustomizeLabelsDialog::toggleNumpadNav);
    connect(okBtn, &QPushButton::clicked, this, &CustomizeLabelsDialog::apply);
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
}

void CustomizeLabelsDialog::toggleNumpadNav(int state) {
    *m_numpadNav = (state == Qt::Checked);
    updateVisibility();
    if (*m_numpadNav && isVisible()) {
        apply();
        m_shouldReopen = true;
        accept();
    }
}

void CustomizeLabelsDialog::updateVisibility() {
    for (int i = 0; i < 10; ++i) {
        bool hidden = (*m_numpadNav) && isNumPadSkip(i);
        m_numLabels[i]->setVisible(!hidden);
        m_edits[i]->setVisible(!hidden);
    }
    adjustSize();
}

void CustomizeLabelsDialog::apply() {
    for (int i = 0; i < 10; ++i) {
        m_labels[i] = m_edits[i]->text();
    }
    emit labelsChanged();
}
