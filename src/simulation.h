#pragma once

#include "graphics/shape.h"
#include <QtCore>
#include "graphics/pointcloud.h"
#include "exporter.h"

class Shader;

struct Fluid {
    std::string name;
    double mass;
    double viscosity;
    double gasConstant;
    float colorI; // interface tension color
    float colorS; // surface tension color
    int numParticles;
    double density;
    double restDensity;
    bool temperatureDependent; // Whether this fluid's density depends on temperature
    double temperatureConstant;
};

struct Particle {
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double pressure;
    double density;
    double restDensity;
    double temperature; // in celsius

    Fluid* fluid;
};

class Simulation
{
public:
    Simulation(QSettings& settings);
    ~Simulation();

    void init();
    void update(double seconds);
    void draw(Shader *shader);

    // Method to update simulation parameters at runtime
    void updateParameters(
        const std::map<QString, float>& paramValues
    );

    // Helper methods for simulation
    double density_S(int i);
    double calculateTemperatureDiffusionStep(int i);
    double calculatePressure(int i, double density, double restDensity);
    Eigen::Vector3d fPressure(int i);
    Eigen::Vector3d fViscosity(int i);
    Eigen::Vector3d fSurfaceTension(int i);
    Eigen::Vector3d fInterfaceTension(int i);
    Eigen::Vector3d checkCollision(Eigen::Vector3d pos);
    void evaluateCollisions(int i);

    void reinitialize();

    // Getter for the number of fluids
    size_t getFluidCount() const { return m_fluids.size(); }
    
    // Add new fluid to the simulation
    void addFluid(const Fluid& fluidProperties);

private:
    Shape m_shape;
    Shape m_box;
    Shape m_collider;
    QSettings& _settings;
    std::vector<Eigen::Vector3d> m_vertices;
    std::vector<Eigen::Vector3i> m_faces;

    std::vector<Particle*> m_particles;
    std::vector<int> m_fluidStartIndices;
    std::vector<std::unique_ptr<Fluid>> m_fluids;
    std::vector<std::unique_ptr<PointCloud>> m_pointclouds;
    float radius;
    float ceiling;
    float coneTop;

    Shape m_ground;
    Exporter m_exporter;

    void initBox();
    void initGround();
    void updateFluidIndices();
    void initFluidsFromSettings();
    void createFluidParticles(Fluid* fluid, const Eigen::Vector3d& startPos, const Eigen::Vector3d& dimensions, double spacing);


    // Global simulation paramters
    Eigen::Vector3d gravity;
    float h; // smoothing length
    float surfaceTensionThreshold;
    float surfaceTensionCoeff;
    float interfaceTensionThreshold;
    float interfaceTensionCoeff;
    float diffusionCoeff;

    // Kernel coefficients
    float Wpoly6Coeff;
    float Wpoly6GradCoeff;
    float Wpoly6LaplacianCoeff;
    float WspikyGradCoeff;
    float WviscosityLaplacianCoeff;

    // Helper methods
    void updateKernelCoefficients();
    double calculateTimeStep();
    void loadFluidParameters();
};

