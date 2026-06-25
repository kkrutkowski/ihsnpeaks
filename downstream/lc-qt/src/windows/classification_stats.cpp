#include "classification_stats.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include <QFrame>

ClassificationStatsDialog::ClassificationStatsDialog(const QString labels[10], const int counts[10], QWidget *parent)
    : QDialog(parent) {

    setWindowTitle("Classification Statistics");
    setMinimumWidth(320);
    setStyleSheet(
        "QDialog { background-color: #242429; color: #e2e8f0; }"
        "QLabel { color: #cbd5e1; }"
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

    QGroupBox *tableBox = new QGroupBox("Counts per Label");
    QGridLayout *grid = new QGridLayout(tableBox);
    grid->setSpacing(6);

    // Header row
    QLabel *hdrNum = new QLabel("#");
    hdrNum->setMinimumWidth(20);
    hdrNum->setAlignment(Qt::AlignCenter);
    QLabel *hdrLabel = new QLabel("Label");
    QLabel *hdrCount = new QLabel("Count");
    hdrCount->setMinimumWidth(50);
    hdrCount->setAlignment(Qt::AlignCenter);

    QFont boldFont;
    boldFont.setBold(true);
    hdrNum->setFont(boldFont);
    hdrLabel->setFont(boldFont);
    hdrCount->setFont(boldFont);

    grid->addWidget(hdrNum, 0, 0);
    grid->addWidget(hdrLabel, 0, 1);
    grid->addWidget(hdrCount, 0, 2);

    for (int i = 0; i < 10; ++i) {
        QLabel *numLabel = new QLabel(QString::number(i));
        numLabel->setMinimumWidth(20);
        numLabel->setAlignment(Qt::AlignCenter);
        grid->addWidget(numLabel, i + 1, 0);

        QLabel *labelText = new QLabel(labels[i]);
        grid->addWidget(labelText, i + 1, 1);

        QLabel *countText = new QLabel(QString::number(counts[i]));
        countText->setMinimumWidth(50);
        countText->setAlignment(Qt::AlignCenter);
        grid->addWidget(countText, i + 1, 2);
    }

    mainLayout->addWidget(tableBox);

    QHBoxLayout *btnLayout = new QHBoxLayout();
    btnLayout->addStretch();
    QPushButton *closeBtn = new QPushButton("Close");
    btnLayout->addWidget(closeBtn);
    mainLayout->addLayout(btnLayout);

    connect(closeBtn, &QPushButton::clicked, this, &QDialog::accept);
}
