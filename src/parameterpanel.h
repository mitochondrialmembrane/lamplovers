#pragma once

#include <QWidget>
#include <QSettings>
#include <QSlider>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <map>
#include <string>

class Simulation;

class ParameterPanel : public QWidget
{
    Q_OBJECT

public:
    ParameterPanel(QSettings& settings, Simulation* simulation, QWidget* parent = nullptr);
    ~ParameterPanel();

private slots:
    void resetToDefaults();
    void updateParameter(int value);

private:
    void createUI();
    void setupConnections();
    void updateSimulationParameters();
    
    // Convert slider value to parameter value (mapping integer to float)
    float sliderToParam(int sliderValue, float minValue, float maxValue);
    // Convert parameter value to slider value (mapping float to integer)
    int paramToSlider(float paramValue, float minValue, float maxValue);

    QSettings& m_settings;
    Simulation* m_simulation;
    
    struct ParameterInfo {
        QString name;
        QString displayName;
        QString tooltip;
        float minValue;
        float maxValue;
        float defaultValue;
        QSlider* slider;
        QLabel* valueLabel;
    };
    
    std::map<QString, ParameterInfo> m_parameters;
    QPushButton* m_resetButton;
};
