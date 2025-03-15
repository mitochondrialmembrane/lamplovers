#include "mainwindow.h"
#include <QHBoxLayout>

MainWindow::MainWindow(QSettings& settings)
{
    glWidget = new GLWidget(settings);

    QHBoxLayout *container = new QHBoxLayout;
    container->addWidget(glWidget);
    this->setLayout(container);
}

MainWindow::~MainWindow()
{
    delete glWidget;
}
