#pragma once

#include <QWidget>
#include <QSettings>
#include <QSlider>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QComboBox>
#include <QTabWidget>
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
    void resetSimulation();
    void updateParameter(int value);
    void fluidSelectionChanged(int index);

private:
    void createUI();
    void setupConnections();
    void updateSimulationParameters();
    void updateFluidParamUi(int fluidIndex);
    
    // Convert slider value to parameter value (mapping integer to float)
    float sliderToParam(int sliderValue, float minValue, float maxValue);
    // Convert parameter value to slider value (mapping float to integer)
    int paramToSlider(float paramValue, float minValue, float maxValue);

    QSettings& m_settings;
    Simulation* m_simulation;
    QTabWidget* m_tabWidget;
    QComboBox* m_fluidSelector;
    
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
    
    // Global parameters (smoothing length, gravity, etc.)
    std::map<QString, ParameterInfo> m_globalParameters;
    
    // Fluid parameters per fluid (density, viscosity, etc.)
    std::vector<std::map<QString, ParameterInfo>> m_fluidParameters;
    
    QPushButton* m_resetButton;
    QPushButton* m_resetSimulationButton;
    QPushButton* m_addFluidButton;
};
