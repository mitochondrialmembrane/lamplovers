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
    
    // Create tab widget for organized parameters
    m_tabWidget = new QTabWidget();
    
    // ===== Create Global Parameters Tab =====
    QWidget* globalParamTab = new QWidget();
    QVBoxLayout* globalLayout = new QVBoxLayout(globalParamTab);
    QGroupBox* globalGroup = new QGroupBox("Simulation Parameters");
    QGridLayout* globalGrid = new QGridLayout(globalGroup);
    
    // Define global parameters
    m_globalParameters["smoothingLength"] = {
        "smoothingLength", 
        "Smoothing Length", 
        "Controls the radius of influence between particles. Affects simulation stability and detail.",
        0.05f, 1.0f,
        m_settings.value("Parameters/smoothingLength", 0.15f).toFloat(), 
        nullptr, nullptr
    };
    
    m_globalParameters["gravity_y"] = {
        "gravity_y", 
        "Gravity Y", 
        "Controls the strength of gravitational force. More negative values create stronger downward pull.",
        -20.0f, 0.0f, 
        m_settings.value("Parameters/gravity_y", -10.0f).toFloat(), 
        nullptr, nullptr
    };
    
    m_globalParameters["surfaceTensionThreshold"] = {
        "surfaceTensionThreshold", 
        "Surface Tension Threshold", 
        "Controls when surface tension is applied. Higher values create more water-droplet behavior.",
        0.1f, 2.0f, 
        m_settings.value("Parameters/surfaceTensionThreshold", 0.5f).toFloat(), 
        nullptr, nullptr
    };
    
    m_globalParameters["surfaceTensionCoeff"] = {
        "surfaceTensionCoeff", 
        "Surface Tension Coefficient", 
        "Controls the strength of surface tension. Higher values create more cohesive fluid surfaces.",
        5.0f, 1000.0f,
        m_settings.value("Parameters/surfaceTensionCoeff", 20.0f).toFloat(), 
        nullptr, nullptr
    };
    
    m_globalParameters["interfaceTensionThreshold"] = {
        "interfaceTensionThreshold",
        "Interface Tension Threshold",
        "Controls when interface tension is applied. Higher values create less separation between fluids",
        0.1f, 2.0f,
        m_settings.value("Parameters/interfaceTensionThreshold", 0.5f).toFloat(),
        nullptr, nullptr
    };
    
    m_globalParameters["interfaceTensionCoeff"] = {
        "interfaceTensionCoeff",
        "Interface Tension Coefficient",
        "Controls the strength of interface tension. Higher values cause more separation between fluids.",
        0.0f, 1000.0f,
        m_settings.value("Parameters/interfaceTensionCoeff", 20.0f).toFloat(),
        nullptr, nullptr
    };
    
    m_globalParameters["diffusionCoeff"] = {
        "diffusionCoeff",
        "Diffusion Coefficient",
        "Controls the strength of diffusion. Higher values cause temperature to diffuse faster.",
        0.0f, 0.01f,
        m_settings.value("Parameters/diffusionCoeff", 0.0001f).toFloat(),
        nullptr, nullptr
    };
    
    // Add global parameters to grid
    int row = 0;
    for (auto& [key, param] : m_globalParameters) {
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
        globalGrid->addWidget(nameLabel, row * 2, 0, 1, 5); // Parameter name spans all columns
        
        // Second row for controls
        globalGrid->addWidget(minLabel, row * 2 + 1, 0);     // Column 0: min value
        globalGrid->addWidget(param.slider, row * 2 + 1, 1); // Column 1: slider
        globalGrid->addWidget(maxLabel, row * 2 + 1, 2);     // Column 2: max value
        globalGrid->addWidget(param.valueLabel, row * 2 + 1, 3); // Column 3: current value 
        
        row++;
    }
    
    // Set column stretches
    globalGrid->setColumnStretch(0, 1);  // Min value - small space
    globalGrid->setColumnStretch(1, 8);  // Slider - most space
    globalGrid->setColumnStretch(2, 1);  // Max value - small space
    globalGrid->setColumnStretch(3, 2);  // Current value - a bit more space for numbers
    
    // Add spacing
    globalGrid->setHorizontalSpacing(10);
    globalGrid->setVerticalSpacing(5);
    
    globalLayout->addWidget(globalGroup);
    
    // ===== Create Fluid Parameters Tab =====
    QWidget* fluidParamTab = new QWidget();
    QVBoxLayout* fluidLayout = new QVBoxLayout(fluidParamTab);
    
    // Add fluid selector
    QHBoxLayout* selectorLayout = new QHBoxLayout();
    QLabel* fluidLabel = new QLabel("Select Fluid:");
    m_fluidSelector = new QComboBox();
    
    // Get number of fluids from settings
    int numFluids = m_settings.value("Parameters/num_fluids", 2).toInt();
    
    for (int i = 0; i < numFluids; i++) {
        QString fluidName = m_settings.value(QString("Fluid%1/name").arg(i+1), 
                                            QString("Fluid %1").arg(i+1)).toString();
        m_fluidSelector->addItem(fluidName);
    }
    
    selectorLayout->addWidget(fluidLabel);
    selectorLayout->addWidget(m_fluidSelector);
    selectorLayout->addStretch();
    
    fluidLayout->addLayout(selectorLayout);
    
    // Create group for fluid parameters
    QGroupBox* fluidGroup = new QGroupBox("Fluid Properties");
    QGridLayout* fluidGrid = new QGridLayout(fluidGroup);
    
    // Initialize fluid parameter maps for each fluid
    m_fluidParameters.resize(numFluids);
    
    // Define standard fluid parameters for each fluid
    for (int fluidIndex = 0; fluidIndex < numFluids; fluidIndex++) {
        QString prefix = QString("Fluid%1/").arg(fluidIndex+1);
        
        m_fluidParameters[fluidIndex]["density"] = {
            QString("fluid%1_density").arg(fluidIndex+1),
            "Density",
            "Controls how densely packed fluid particles are. Higher values create heavier fluid.",
            200.0f, 4000.0f,
            m_settings.value(prefix + "density", 1000.0f).toFloat(),
            nullptr, nullptr
        };
        
        m_fluidParameters[fluidIndex]["viscosity"] = {
            QString("fluid%1_viscosity").arg(fluidIndex+1), 
            "Viscosity",
            "Controls how thick/sticky the fluid is. Higher values create more honey-like fluid.",
            1.0f, 1000.0f,
            m_settings.value(prefix + "viscosity", 40.0f).toFloat(),
            nullptr, nullptr
        };
        
        m_fluidParameters[fluidIndex]["idealGasConstant"] = {
            QString("fluid%1_idealGasConstant").arg(fluidIndex+1),
            "Gas Constant",
            "Controls pressure response. Higher values create more repulsive force between particles.",
            10.0f, 100.0f, 
            m_settings.value(prefix + "idealGasConstant", 40.0f).toFloat(),
            nullptr, nullptr
        };
        
        // Add color parameters if needed
        m_fluidParameters[fluidIndex]["colorI"] = {
            QString("fluid%1_colorI").arg(fluidIndex+1),
            "Interface Color",
            "Controls interface tension between fluids. Range from -1 to 1.",
            -1.0f, 1.0f,
            m_settings.value(prefix + "colorI", fluidIndex == 0 ? -0.5f : 0.5f).toFloat(),
            nullptr, nullptr
        };
    }
    
    // Create UI elements for first fluid parameters
    row = 0;
    for (auto& [key, param] : m_fluidParameters[0]) {
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
        fluidGrid->addWidget(nameLabel, row * 2, 0, 1, 5); // Parameter name spans all columns
        
        // Second row for controls
        fluidGrid->addWidget(minLabel, row * 2 + 1, 0);     // Column 0: min value
        fluidGrid->addWidget(param.slider, row * 2 + 1, 1); // Column 1: slider
        fluidGrid->addWidget(maxLabel, row * 2 + 1, 2);     // Column 2: max value
        fluidGrid->addWidget(param.valueLabel, row * 2 + 1, 3); // Column 3: current value 
        
        row++;
    }
    
    // Set column stretches
    fluidGrid->setColumnStretch(0, 1);  // Min value - small space
    fluidGrid->setColumnStretch(1, 8);  // Slider - most space
    fluidGrid->setColumnStretch(2, 1);  // Max value - small space
    fluidGrid->setColumnStretch(3, 2);  // Current value - a bit more space for numbers
    
    // Add spacing
    fluidGrid->setHorizontalSpacing(10);
    fluidGrid->setVerticalSpacing(5);
    
    fluidLayout->addWidget(fluidGroup);
    fluidLayout->addStretch();
    
    // Add tabs to the tab widget
    m_tabWidget->addTab(globalParamTab, "Global Parameters");
    m_tabWidget->addTab(fluidParamTab, "Fluid Parameters");
    
    // Create buttons
    m_resetButton = new QPushButton("Reset to Defaults");
    m_resetSimulationButton = new QPushButton("Reset Simulation");
    
    // Create a button layout
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->addWidget(m_resetButton);
    buttonLayout->addWidget(m_resetSimulationButton);
    
    // Add everything to main layout
    mainLayout->addWidget(m_tabWidget);
    mainLayout->addLayout(buttonLayout);
    
    // Set layout margins
    mainLayout->setContentsMargins(5, 5, 5, 5);
    
    setMinimumWidth(350);
    setMaximumWidth(450);
}

void ParameterPanel::setupConnections()
{
    // Connect reset buttons
    connect(m_resetButton, &QPushButton::clicked, this, &ParameterPanel::resetToDefaults);
    connect(m_resetSimulationButton, &QPushButton::clicked, this, &ParameterPanel::resetSimulation);
    
    // Connect fluid selector
    connect(m_fluidSelector, QOverload<int>::of(&QComboBox::currentIndexChanged), 
            this, &ParameterPanel::fluidSelectionChanged);
    
    // Connect global parameter sliders
    for (auto& [key, param] : m_globalParameters) {
        connect(param.slider, &QSlider::valueChanged, this, &ParameterPanel::updateParameter);
    }
    
    // Connect fluid parameter sliders for the first fluid (others are connected as needed)
    for (auto& [key, param] : m_fluidParameters[0]) {
        connect(param.slider, &QSlider::valueChanged, this, &ParameterPanel::updateParameter);
    }
}

void ParameterPanel::resetToDefaults()
{
    // Reset all global sliders to default values
    for (auto& [key, param] : m_globalParameters) {
        param.slider->setValue(paramToSlider(param.defaultValue, param.minValue, param.maxValue));
    }
    
    // Reset fluid sliders for all fluids
    for (size_t fluidIdx = 0; fluidIdx < m_fluidParameters.size(); fluidIdx++) {
        for (auto& [key, param] : m_fluidParameters[fluidIdx]) {
            param.slider->setValue(paramToSlider(param.defaultValue, param.minValue, param.maxValue));
        }
    }
    
    // Update simulation with default parameters
    updateSimulationParameters();
}

void ParameterPanel::updateParameter(int value) {
    QObject* sender = QObject::sender();
    if (!sender)
        return;
        
    QString name = sender->objectName();
    
    // Find parameter and update its value label (in global parameters)
    bool found = false;
    for (auto& [key, param] : m_globalParameters) {
        if (param.name == name && param.slider) {  // Add null check
            float paramValue = sliderToParam(value, param.minValue, param.maxValue);
            param.valueLabel->setText(QString::number(paramValue, 'f', 2));
            found = true;
            break;
        }
    }
    
    // If not found in global parameters, check fluid parameters
    if (!found) {
        for (size_t fluidIdx = 0; fluidIdx < m_fluidParameters.size(); fluidIdx++) {
            for (auto& [key, param] : m_fluidParameters[fluidIdx]) {
                if (param.name == name && param.slider) {  // Add null check
                    float paramValue = sliderToParam(value, param.minValue, param.maxValue);
                    if (param.valueLabel) {  // Add null check 
                        param.valueLabel->setText(QString::number(paramValue, 'f', 2));
                    }
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
    }
    
    // Update all simulation parameters
    updateSimulationParameters();
}

void ParameterPanel::updateSimulationParameters() {
    // Create a map of parameter names to values
    std::map<QString, float> paramValues;
    
    // Add global parameters
    for (const auto& [key, param] : m_globalParameters) {
        if (param.slider) {  // Add null check
            float value = sliderToParam(param.slider->value(), param.minValue, param.maxValue);
            paramValues[param.name] = value;
        }
    }
    
    // Add fluid parameters for all fluids
    for (size_t fluidIdx = 0; fluidIdx < m_fluidParameters.size(); fluidIdx++) {
        for (const auto& [key, param] : m_fluidParameters[fluidIdx]) {
            if (param.slider) {  // Add null check
                float value = sliderToParam(param.slider->value(), param.minValue, param.maxValue);
                paramValues[param.name] = value;
            }
        }
    }
    
    // Update the simulation with new parameter values
    if (m_simulation) {  // Add null check
        m_simulation->updateParameters(paramValues);
    }
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
    try {
        // First reset parameters to a safe state
        resetToDefaults();
        
        // Then tell the simulation to reinitialize with a try/catch block
        if (m_simulation) {
            m_simulation->reinitialize();
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error resetting simulation: " << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown error resetting simulation" << std::endl;
    }
}

void ParameterPanel::fluidSelectionChanged(int index)
{
    // Update the UI to show the parameters for the selected fluid
    updateFluidParamUi(index);
}

void ParameterPanel::updateFluidParamUi(int fluidIndex)
{
    // Make sure the index is valid
    if (fluidIndex < 0 || fluidIndex >= static_cast<int>(m_fluidParameters.size())) {
        return;
    }
    
    // Disconnect all existing fluid parameter sliders
    for (size_t i = 0; i < m_fluidParameters.size(); i++) {
        for (auto& [key, param] : m_fluidParameters[i]) {
            if (param.slider) {
                disconnect(param.slider, &QSlider::valueChanged, this, &ParameterPanel::updateParameter);
            }
        }
    }
    
    // Find the fluid parameters group
    QGroupBox* fluidGroup = nullptr;
    QGridLayout* fluidGrid = nullptr;
    
    if (m_tabWidget->count() > 1) {
        QWidget* fluidTab = m_tabWidget->widget(1); // Fluid parameters tab
        if (fluidTab) {
            QLayout* layout = fluidTab->layout();
            if (layout) {
                // Find the fluid parameters group box
                for (int i = 0; i < layout->count(); i++) {
                    QLayoutItem* item = layout->itemAt(i);
                    if (item && item->widget() && item->widget()->inherits("QGroupBox")) {
                        fluidGroup = qobject_cast<QGroupBox*>(item->widget());
                        if (fluidGroup) {
                            fluidGrid = qobject_cast<QGridLayout*>(fluidGroup->layout());
                            break;
                        }
                    }
                }
            }
        }
    }
    
    if (!fluidGroup || !fluidGrid) {
        return;
    }
    
    // Clear the existing fluid parameter widgets
    QLayoutItem* child;
    while ((child = fluidGrid->takeAt(0)) != nullptr) {
        if (child->widget()) {
            child->widget()->hide();
            delete child->widget();
        }
        delete child;
    }
    
    // Add the selected fluid's parameters
    int row = 0;
    for (auto& [key, param] : m_fluidParameters[fluidIndex]) {
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
        fluidGrid->addWidget(nameLabel, row * 2, 0, 1, 5); // Parameter name spans all columns
        
        // Second row for controls
        fluidGrid->addWidget(minLabel, row * 2 + 1, 0);     // Column 0: min value
        fluidGrid->addWidget(param.slider, row * 2 + 1, 1); // Column 1: slider
        fluidGrid->addWidget(maxLabel, row * 2 + 1, 2);     // Column 2: max value
        fluidGrid->addWidget(param.valueLabel, row * 2 + 1, 3); // Column 3: current value 
        
        // Connect slider
        connect(param.slider, &QSlider::valueChanged, this, &ParameterPanel::updateParameter);
        
        row++;
    }
    
    // Set column stretches
    fluidGrid->setColumnStretch(0, 1);  // Min value - small space
    fluidGrid->setColumnStretch(1, 8);  // Slider - most space
    fluidGrid->setColumnStretch(2, 1);  // Max value - small space
    fluidGrid->setColumnStretch(3, 2);  // Current value - a bit more space for numbers
    
    // Add spacing
    fluidGrid->setHorizontalSpacing(10);
    fluidGrid->setVerticalSpacing(5);
}
