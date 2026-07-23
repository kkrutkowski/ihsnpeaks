#pragma once

#include <QDialog>
#include <QDoubleSpinBox>

struct PeriodSearchSettings;

/*
 * PeriodScrollDialog — edits the period scroll rate: how fast the "<" / ">"
 * arrow buttons slide the pivot frequency while held, expressed in 1/DeltaT
 * per second (DeltaT = time from the first to the last measurement).
 */
class PeriodScrollDialog : public QDialog {
    Q_OBJECT
public:
    explicit PeriodScrollDialog(PeriodSearchSettings *settings, QWidget *parent = nullptr);

private slots:
    void apply();

private:
    PeriodSearchSettings *m_settings;
    QDoubleSpinBox *m_rate;
};
