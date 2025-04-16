#pragma once

#include <QMainWindow>
#include "glwidget.h"
#include "parameterpanel.h"
#include <QtCore>

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    MainWindow(QSettings& settings);
    ~MainWindow();

private:
    GLWidget *glWidget;
    ParameterPanel *parameterPanel;
};
