#pragma once

#include <QFrame>
#include <QPainter>
#include <QVector>
#include <QString>
#include <cmath>

extern "C" {
#include "lc_readout.h"
}

class LightCurvePlotWidget : public QFrame {
    Q_OBJECT
public:
    explicit LightCurvePlotWidget(const QString &title, QWidget *parent = nullptr);
    ~LightCurvePlotWidget() override;

    void setData(const lc_data_t *data);
    void clearData();
    bool hasData() const { return m_n > 0; }

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    QString m_title;
    QVector<double> m_x;
    QVector<float> m_y;
    unsigned int m_n = 0;
};
