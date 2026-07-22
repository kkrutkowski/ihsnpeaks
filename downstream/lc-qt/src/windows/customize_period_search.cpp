#include "customize_period_search.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QLocale>
#include <QDoubleValidator>

CustomizePeriodSearchDialog::CustomizePeriodSearchDialog(PeriodSearchSettings *settings, QWidget *parent)
    : QDialog(parent), m_settings(settings) {

    setWindowTitle("Customize Period Search");
    setMinimumWidth(360);
    setStyleSheet(
        "QDialog { background-color: #242429; color: #e2e8f0; }"
        "QLabel { color: #cbd5e1; }"
        "QSpinBox, QDoubleSpinBox, QComboBox, QLineEdit { background-color: #17171b; border: 1px solid #32323b; "
        "  border-radius: 4px; color: #ffffff; padding: 3px; }"
        "QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QLineEdit:focus { border: 1px solid #3b82f6; }"
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

    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    mainLayout->setContentsMargins(12, 12, 12, 12);
    mainLayout->setSpacing(8);

    QGroupBox *box = new QGroupBox("Period Search Parameters");
    QGridLayout *grid = new QGridLayout(box);
    grid->setSpacing(6);

    int row = 0;

    m_terms = new QSpinBox();
    m_terms->setRange(1, 32);
    m_terms->setValue(m_settings->nterms);
    grid->addWidget(new QLabel("Number of harmonics"), row, 0);
    grid->addWidget(m_terms, row++, 1);

    m_oversampling = new QDoubleSpinBox();
    m_oversampling->setLocale(QLocale::c());
    m_oversampling->setRange(0.1, 1000.0);
    m_oversampling->setDecimals(2);
    m_oversampling->setSingleStep(0.5);
    m_oversampling->setValue(m_settings->oversampling);
    grid->addWidget(new QLabel("Oversampling"), row, 0);
    grid->addWidget(m_oversampling, row++, 1);

    m_fmin = new QLineEdit();
    m_fmin->setPlaceholderText("auto");
    m_fmin->setValidator(new QDoubleValidator(0.0, 1e6, 6, m_fmin));
    if (m_settings->fmin > 0.0)
        m_fmin->setText(QString::number(m_settings->fmin, 'g', 10));
    grid->addWidget(new QLabel("Min frequency"), row, 0);
    grid->addWidget(m_fmin, row++, 1);

    m_fmax = new QDoubleSpinBox();
    m_fmax->setLocale(QLocale::c());
    m_fmax->setRange(0.0001, 1e6);
    m_fmax->setDecimals(4);
    m_fmax->setValue(m_settings->fmax);
    grid->addWidget(new QLabel("Max frequency"), row, 0);
    grid->addWidget(m_fmax, row++, 1);

    m_zoom = new QDoubleSpinBox();
    m_zoom->setLocale(QLocale::c());
    m_zoom->setRange(1.0, 1000.0);
    m_zoom->setDecimals(2);
    m_zoom->setValue(m_settings->zoomFactor);
    grid->addWidget(new QLabel("Zooming factor"), row, 0);
    grid->addWidget(m_zoom, row++, 1);

    m_radius = new QDoubleSpinBox();
    m_radius->setLocale(QLocale::c());
    m_radius->setRange(0.0, 100.0);
    m_radius->setDecimals(4);
    m_radius->setValue(m_settings->searchRadius);
    grid->addWidget(new QLabel("Search radius"), row, 0);
    grid->addWidget(m_radius, row++, 1);

    m_pswf = new QComboBox();
    m_pswf->addItem("4/3", 43);
    m_pswf->addItem("2/1", 21);
    m_pswf->setCurrentIndex(m_settings->pswf == 21 ? 1 : 0);
    grid->addWidget(new QLabel("Upsampling factor"), row, 0);
    grid->addWidget(m_pswf, row++, 1);

    m_oversmoothing = new QDoubleSpinBox();
    m_oversmoothing->setLocale(QLocale::c());
    m_oversmoothing->setRange(0.0, 1.0);
    m_oversmoothing->setDecimals(4);
    m_oversmoothing->setSingleStep(0.05);
    m_oversmoothing->setValue(m_settings->oversmoothing);
    grid->addWidget(new QLabel("Oversmoothing factor"), row, 0);
    grid->addWidget(m_oversmoothing, row++, 1);

    m_nbins = new QSpinBox();
    m_nbins->setRange(1, 10000);
    m_nbins->setValue(m_settings->nbins);
    grid->addWidget(new QLabel("Number of bins (BLS)"), row, 0);
    grid->addWidget(m_nbins, row++, 1);

    m_autoCenter = new QCheckBox("Automatically center on highest peak");
    m_autoCenter->setChecked(m_settings->autoCenter);
    grid->addWidget(m_autoCenter, row++, 0, 1, 2);

    grid->setColumnStretch(1, 1);
    mainLayout->addWidget(box);

    QHBoxLayout *btnLayout = new QHBoxLayout();
    btnLayout->addStretch();
    QPushButton *okBtn = new QPushButton("OK");
    QPushButton *cancelBtn = new QPushButton("Cancel");
    btnLayout->addWidget(okBtn);
    btnLayout->addWidget(cancelBtn);
    mainLayout->addLayout(btnLayout);

    connect(okBtn, &QPushButton::clicked, this, &CustomizePeriodSearchDialog::apply);
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
}

void CustomizePeriodSearchDialog::apply() {
    m_settings->nterms = m_terms->value();
    m_settings->oversampling = m_oversampling->value();
    m_settings->fmin = m_fmin->text().isEmpty() ? 0.0 : m_fmin->text().toDouble();
    m_settings->fmax = m_fmax->value();
    m_settings->zoomFactor = m_zoom->value();
    m_settings->searchRadius = m_radius->value();
    m_settings->pswf = m_pswf->currentData().toInt();
    m_settings->oversmoothing = m_oversmoothing->value();
    m_settings->nbins = m_nbins->value();
    m_settings->autoCenter = m_autoCenter->isChecked();
}
