#include <QApplication>
#include "windows/main_window.h"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    app.setStyle("Fusion");

    MainWindow window;
    window.show();

    if (argc > 1) {
        window.loadLightCurve(QString::fromLocal8Bit(argv[1]));
    }

    return app.exec();
}
