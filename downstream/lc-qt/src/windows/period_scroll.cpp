#include "period_scroll.h"
#include "customize_period_search.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QLocale>

PeriodScrollDialog::PeriodScrollDialog(PeriodSearchSettings *settings, QWidget *parent)
    : QDialog(parent), m_settings(settings) {

    setWindowTitle("Period Scroll");
    setMinimumWidth(320);
    setStyleSheet(
        "QDialog { background-color: #242429; color: #e2e8f0; }"
        "QLabel { color: #cbd5e1; }"
        "QDoubleSpinBox { background-color: #17171b; border: 1px solid #32323b; "
        "  border-radius: 4px; color: #ffffff; padding: 3px; }"
        "QDoubleSpinBox:focus { border: 1px solid #3b82f6; }"
        "QPushButton { background-color: #303038; border: 1px solid #3d3d47; border-radius: 4px; "
        "  color: #f1f5f9; padding: 4px 12px; font-weight: bold; }"
        "QPushButton:hover { background-color: #3b82f6; border-color: #60a5fa; }"
        "QPushButton:pressed { background-color: #1d4ed8; }"
        "QGroupBox { border: 1px solid #32323b; border-radius: 6px; margin-top: 10px; padding: 6px; "
        "  font-weight: bold; color: #60a5fa; }"
        "QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; "
        "  background-color: #242429; }"
    );

    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    mainLayout->setContentsMargins(12, 12, 12, 12);
    mainLayout->setSpacing(8);

    QGroupBox *box = new QGroupBox("Period Scroll");
    QGridLayout *grid = new QGridLayout(box);
    grid->setSpacing(6);

    m_rate = new QDoubleSpinBox();
    m_rate->setLocale(QLocale::c());
    m_rate->setRange(0.01, 1000.0);
    m_rate->setDecimals(2);
    m_rate->setSingleStep(0.1);
    m_rate->setValue(m_settings->scrollRate);
    grid->addWidget(new QLabel("Scroll rate (1/\u0394T per second)"), 0, 0);
    grid->addWidget(m_rate, 0, 1);

    grid->setColumnStretch(1, 1);
    mainLayout->addWidget(box);

    QHBoxLayout *btnLayout = new QHBoxLayout();
    btnLayout->addStretch();
    QPushButton *okBtn = new QPushButton("OK");
    QPushButton *cancelBtn = new QPushButton("Cancel");
    btnLayout->addWidget(okBtn);
    btnLayout->addWidget(cancelBtn);
    mainLayout->addLayout(btnLayout);

    connect(okBtn, &QPushButton::clicked, this, &PeriodScrollDialog::apply);
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
}

void PeriodScrollDialog::apply() {
    m_settings->scrollRate = m_rate->value();
}
