#include "mainwindow.h"
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QSplitter>

MainWindow::MainWindow(QSettings& settings)
{
    // Create GLWidget
    glWidget = new GLWidget(settings);
    
    // Create parameter panel
    parameterPanel = new ParameterPanel(settings, &glWidget->getSimulation());
    
    // Create splitter for resizable panels
    QSplitter *splitter = new QSplitter(Qt::Horizontal);
    splitter->addWidget(parameterPanel);
    splitter->addWidget(glWidget);
    
    // Set initial sizes
    splitter->setSizes(QList<int>() << 300 << 600);
    
    // Add to main layout
    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(splitter);
    layout->setContentsMargins(0, 0, 0, 0);
    setLayout(layout);
}

MainWindow::~MainWindow()
{
}
