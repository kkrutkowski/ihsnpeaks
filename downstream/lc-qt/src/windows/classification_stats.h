#pragma once

#include <QDialog>
#include <QString>

class ClassificationStatsDialog : public QDialog {
    Q_OBJECT
public:
    explicit ClassificationStatsDialog(const QString labels[10], const int counts[10], QWidget *parent = nullptr);
};
