#include "parameterpanel.h"
#include "simulation.h"

#include <QFormLayout>
#include <QDoubleSpinBox>
#include <QDebug>

ParameterPanel::ParameterPanel(QSettings& settings, Simulation* simulation, QWidget* parent)
    : QWidget(parent)
    , m_settings(settings)
    , m_simulation(simulation)
{
    createUI();
    setupConnections();
}

ParameterPanel::~ParameterPanel()
{
}

void ParameterPanel::createUI()
{
    // Create main layout
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QGroupBox* paramGroup = new QGroupBox("Simulation Parameters");
    QFormLayout* formLayout = new QFormLayout(paramGroup);
    
    // Define parameters to create sliders for
    m_parameters["fluid1_density"] = {
        "fluid1_density", "Fluid Density", 500.0f, 2000.0f, 
        m_settings.value("Parameters/fluid1_density", 1000.0f).toFloat(), nullptr, nullptr
    };
    
    m_parameters["fluid1_viscosity"] = {
        "fluid1_viscosity", "Fluid Viscosity", 1.0f, 20.0f, 
        m_settings.value("Parameters/fluid1_viscosity", 9.0f).toFloat(), nullptr, nullptr
    };
    
    m_parameters["gravity_y"] = {
        "gravity_y", "Gravity Y", -20.0f, 0.0f, 
        m_settings.value("Parameters/gravity_y", -10.0f).toFloat(), nullptr, nullptr
    };
    
    m_parameters["smoothingLength"] = {
        "smoothingLength", "Smoothing Length", 0.05f, 0.3f, 
        m_settings.value("Parameters/smoothingLength", 0.15f).toFloat(), nullptr, nullptr
    };
    
    m_parameters["idealGasConstant"] = {
        "idealGasConstant", "Gas Constant", 10.0f, 100.0f, 
        m_settings.value("Parameters/idealGasConstant", 40.0f).toFloat(), nullptr, nullptr
    };
    
    m_parameters["surfaceTensionThreshold"] = {
        "surfaceTensionThreshold", "Surface Tension Threshold", 0.1f, 2.0f, 
        m_settings.value("Parameters/surfaceTensionThreshold", 0.5f).toFloat(), nullptr, nullptr
    };
    
    m_parameters["surfaceTensionCoeff"] = {
        "surfaceTensionCoeff", "Surface Tension Coefficient", 5.0f, 50.0f, 
        m_settings.value("Parameters/surfaceTensionCoeff", 20.0f).toFloat(), nullptr, nullptr
    };
    
    // Create sliders and labels for each parameter
    for (auto& [key, param] : m_parameters) {
        QHBoxLayout* sliderLayout = new QHBoxLayout();
        
        // Create and configure slider
        param.slider = new QSlider(Qt::Horizontal);
        param.slider->setRange(0, 1000);
        param.slider->setValue(paramToSlider(param.defaultValue, param.minValue, param.maxValue));
        param.slider->setObjectName(param.name);
        param.slider->setMinimumWidth(200);
        
        // Create value label
        param.valueLabel = new QLabel(QString::number(param.defaultValue, 'f', 2));
        param.valueLabel->setMinimumWidth(50);
        param.valueLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        
        // Add to layout
        sliderLayout->addWidget(param.slider);
        sliderLayout->addWidget(param.valueLabel);
        
        formLayout->addRow(param.displayName + ":", sliderLayout);
    }
    
    // Add reset button
    m_resetButton = new QPushButton("Reset to Defaults");
    mainLayout->addWidget(paramGroup);
    mainLayout->addWidget(m_resetButton);
    
    // Set layout margins
    mainLayout->setContentsMargins(5, 5, 5, 5);
    
    setMinimumWidth(300);
    setMaximumWidth(400);
}

void ParameterPanel::setupConnections()
{
    // Connect reset button
    connect(m_resetButton, &QPushButton::clicked, this, &ParameterPanel::resetToDefaults);
    
    // Connect sliders
    for (auto& [key, param] : m_parameters) {
        connect(param.slider, &QSlider::valueChanged, this, &ParameterPanel::updateParameter);
    }
}

void ParameterPanel::resetToDefaults()
{
    // Reset all sliders to default values
    for (auto& [key, param] : m_parameters) {
        param.slider->setValue(paramToSlider(param.defaultValue, param.minValue, param.maxValue));
    }
    
    // Update simulation with default parameters
    updateSimulationParameters();
}

void ParameterPanel::updateParameter(int value)
{
    QObject* sender = QObject::sender();
    if (!sender)
        return;
        
    QString name = sender->objectName();
    
    // Find parameter and update its value label
    for (auto& [key, param] : m_parameters) {
        if (param.name == name) {
            float paramValue = sliderToParam(value, param.minValue, param.maxValue);
            param.valueLabel->setText(QString::number(paramValue, 'f', 2));
            break;
        }
    }
    
    // Update all simulation parameters
    updateSimulationParameters();
}

void ParameterPanel::updateSimulationParameters()
{
    // Create a map of parameter names to values
    std::map<QString, float> paramValues;
    for (const auto& [key, param] : m_parameters) {
        float value = sliderToParam(param.slider->value(), param.minValue, param.maxValue);
        paramValues[param.name] = value;
    }
    
    // Update the simulation with new parameter values
    m_simulation->updateParameters(
        paramValues["fluid1_density"],
        paramValues["fluid1_viscosity"],
        Eigen::Vector3d(0, paramValues["gravity_y"], 0),
        paramValues["smoothingLength"],
        paramValues["idealGasConstant"],
        paramValues["surfaceTensionThreshold"],
        paramValues["surfaceTensionCoeff"]
    );
}

float ParameterPanel::sliderToParam(int sliderValue, float minValue, float maxValue)
{
    return minValue + (sliderValue / 1000.0f) * (maxValue - minValue);
}

int ParameterPanel::paramToSlider(float paramValue, float minValue, float maxValue)
{
    return static_cast<int>((paramValue - minValue) / (maxValue - minValue) * 1000.0f);
}
