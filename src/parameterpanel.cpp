#include "parameterpanel.h"
#include "simulation.h"

#include <QFormLayout>
#include <QDoubleSpinBox>
#include <QDebug>
#include <iostream>

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
    // QFormLayout* formLayout = new QFormLayout(paramGroup);
    
    // Define parameters to create sliders for
    m_parameters["fluid1_density"] = {
        "fluid1_density",
        "Fluid 1 Density",
        "Controls how densely packed fluid particles are. Higher values create heavier fluid.",
        200.0f, 4000.0f,
        m_settings.value("Parameters/fluid1_density", 1000.0f).toFloat(), 
        nullptr, nullptr
    };

    m_parameters["fluid2_density"] = {
        "fluid2_density",
        "Fluid 2 Density",
        "Controls how densely packed fluid particles are. Higher values create heavier fluid.",
        200.0f, 4000.0f,
        m_settings.value("Parameters/fluid2_density", 2000.0f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["fluid1_viscosity"] = {
        "fluid1_viscosity", 
        "Fluid 1 Viscosity",
        "Controls how thick/sticky the fluid is. Higher values create more honey-like fluid.",
        1.0f, 1000.0f,
        m_settings.value("Parameters/fluid1_viscosity", 40.0f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["fluid2_viscosity"] = {
        "fluid2_viscosity",
        "Fluid 2 Viscosity",
        "Controls how thick/sticky the fluid is. Higher values create more honey-like fluid.",
        1.0f, 1000.0f,
        m_settings.value("Parameters/fluid2_viscosity", 40.0f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["gravity_y"] = {
        "gravity_y", 
        "Gravity Y", 
        "Controls the strength of gravitational force. More negative values create stronger downward pull.",
        -20.0f, 0.0f, 
        m_settings.value("Parameters/gravity_y", -10.0f).toFloat(), 
        nullptr, nullptr
    };

    m_parameters["smoothingLength"] = {
        "smoothingLength", 
        "Smoothing Length", 
        "Controls the radius of influence between particles. Affects simulation stability and detail.",
        0.05f, 1.0f,
        m_settings.value("Parameters/smoothingLength", 0.15f).toFloat(), 
        nullptr, nullptr
    };

    m_parameters["fluid1_idealGasConstant"] = {
        "fluid1_idealGasConstant",
        "Fluid 1 Gas Constant",
        "Controls pressure response. Higher values create more repulsive force between particles.",
        10.0f, 100.0f, 
        m_settings.value("Parameters/fluid1_idealGasConstant", 40.0f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["fluid2_idealGasConstant"] = {
        "fluid2_idealGasConstant",
        "Fluid 2 Gas Constant",
        "Controls pressure response. Higher values create more repulsive force between particles.",
        10.0f, 100.0f,
        m_settings.value("Parameters/fluid2_idealGasConstant", 40.0f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["surfaceTensionThreshold"] = {
        "surfaceTensionThreshold", 
        "Surface Tension Threshold", 
        "Controls when surface tension is applied. Higher values create more water-droplet behavior.",
        0.1f, 2.0f, 
        m_settings.value("Parameters/surfaceTensionThreshold", 0.5f).toFloat(), 
        nullptr, nullptr
    };

    m_parameters["surfaceTensionCoeff"] = {
        "surfaceTensionCoeff", 
        "Surface Tension Coefficient", 
        "Controls the strength of surface tension. Higher values create more cohesive fluid surfaces.",
        5.0f, 1000.0f,
        m_settings.value("Parameters/surfaceTensionCoeff", 20.0f).toFloat(), 
        nullptr, nullptr
    };

    m_parameters["interfaceTensionThreshold"] = {
        "interfaceTensionThreshold",
        "Interface Tension Threshold",
        "Controls when interface tension is applied. Higher values create less separation between fluids",
        0.1f, 2.0f,
        m_settings.value("Parameters/interfaceTensionThreshold", 0.5f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["interfaceTensionCoeff"] = {
        "interfaceTensionCoeff",
        "Interface Tension Coefficient",
        "Controls the strength of interface tension. Higher values cause more separation between fluids.",
        0.0f, 1000.0f,
        m_settings.value("Parameters/interfaceTensionCoeff", 20.0f).toFloat(),
        nullptr, nullptr
    };

    m_parameters["diffusionCoeff"] = {
        "diffusionCoeff",
        "Diffusion Coefficient",
        "Controls the strength of diffusion. Higher values cause temperature to diffuse faster.",
        0.0f, 0.01f,
        m_settings.value("Parameters/interfaceTensionCoeff", 0.0001f).toFloat(),
        nullptr, nullptr
    };



    // Add reset button
    QGridLayout* gridLayout = new QGridLayout(paramGroup);
    paramGroup->setLayout(gridLayout);
    int row = 0;

    // Create sliders and labels for each parameter
    for (auto& [key, param] : m_parameters) {
        // Create parameter name label
        QLabel* nameLabel = new QLabel(param.displayName + ":");
        nameLabel->setFont(QFont("", -1, QFont::Bold)); // Make it bold for emphasis
        
        // Create min/max labels
        QLabel* minLabel = new QLabel(QString::number(param.minValue, 'f', 1));
        QLabel* maxLabel = new QLabel(QString::number(param.maxValue, 'f', 1));
        
        // Create slider
        param.slider = new QSlider(Qt::Horizontal);
        param.slider->setRange(0, 1000);
        param.slider->setValue(paramToSlider(param.defaultValue, param.minValue, param.maxValue));
        param.slider->setObjectName(param.name);
        param.slider->setToolTip(param.tooltip);
        param.slider->setMinimumWidth(150);
        
        // Create value label
        param.valueLabel = new QLabel(QString::number(param.defaultValue, 'f', 2));
        param.valueLabel->setMinimumWidth(50);
        param.valueLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        
        // Add to grid layout - use 2 rows per parameter
        gridLayout->addWidget(nameLabel, row * 2, 0, 1, 5); // Parameter name spans all columns

        // Second row for controls - use explicit columns without spanning
        gridLayout->addWidget(minLabel, row * 2 + 1, 0);     // Column 0: min value
        gridLayout->addWidget(param.slider, row * 2 + 1, 1); // Column 1: slider
        gridLayout->addWidget(maxLabel, row * 2 + 1, 2);     // Column 2: max value
        gridLayout->addWidget(param.valueLabel, row * 2 + 1, 3); // Column 3: current value 
        
        row++;
    }

    // Set column stretches - for the new column arrangement
    gridLayout->setColumnStretch(0, 1);  // Min value - small space
    gridLayout->setColumnStretch(1, 8);  // Slider - most space
    gridLayout->setColumnStretch(2, 1);  // Max value - small space
    gridLayout->setColumnStretch(3, 2);  // Current value - a bit more space for numbers
    // Add spacing
    gridLayout->setHorizontalSpacing(10);
    gridLayout->setVerticalSpacing(5);
    
    m_resetButton = new QPushButton("Reset to Defaults");
    m_resetSimulationButton = new QPushButton("Reset Simulation");

    // Create a button layout
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->addWidget(m_resetButton);
    buttonLayout->addWidget(m_resetSimulationButton);

    mainLayout->addWidget(paramGroup);
    mainLayout->addLayout(buttonLayout);
    
    // Set layout margins
    mainLayout->setContentsMargins(5, 5, 5, 5);
    
    setMinimumWidth(300);
    setMaximumWidth(400);
}

void ParameterPanel::setupConnections()
{
    // Connect reset button
    connect(m_resetButton, &QPushButton::clicked, this, &ParameterPanel::resetToDefaults);
    connect(m_resetSimulationButton, &QPushButton::clicked, this, &ParameterPanel::resetSimulation);
    
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
        std::cout << param.name.toStdString() << std::endl;
        std::cout << value << std::endl;
    }
    
    // Update the simulation with new parameter values
    m_simulation->updateParameters(
        paramValues["fluid1_density"],
        paramValues["fluid2_density"],
        paramValues["fluid1_viscosity"],
        paramValues["fluid2_viscosity"],
        Eigen::Vector3d(0, paramValues["gravity_y"], 0),
        paramValues["smoothingLength"],
        paramValues["fluid1_idealGasConstant"],
        paramValues["fluid2_idealGasConstant"],
        paramValues["surfaceTensionThreshold"],
        paramValues["surfaceTensionCoeff"],
        paramValues["interfaceTensionThreshold"],
        paramValues["interfaceTensionCoeff"],
        paramValues["diffusionCoeff"]
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

void ParameterPanel::resetSimulation()
{
    // First reset parameters to defaults
    // resetToDefaults();
    
    // Then tell the simulation to reinitialize
    m_simulation->reinitialize();
}
