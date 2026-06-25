#pragma once

#include <QDialog>
#include <QLineEdit>
#include <QLabel>
#include <QCheckBox>
#include <QGroupBox>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QString>

class CustomizeLabelsDialog : public QDialog {
    Q_OBJECT
public:
    explicit CustomizeLabelsDialog(QString labels[10], bool *numpadNav, QWidget *parent = nullptr);
    bool shouldReopen() const { return m_shouldReopen; }

signals:
    void labelsChanged();

private slots:
    void apply();
    void toggleNumpadNav(int state);

private:
    QString *m_labels;
    bool *m_numpadNav;
    QLineEdit *m_edits[10];
    QCheckBox *m_numpadCheck;
    QLabel *m_numLabels[10];
    QGroupBox *m_gridBox;
    QGridLayout *m_grid;
    QVBoxLayout *m_mainLayout;
    bool m_shouldReopen;

    void updateVisibility();
};
