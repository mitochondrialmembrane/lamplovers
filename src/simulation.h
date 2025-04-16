#pragma once

#include "graphics/shape.h"
#include <QtCore>
#include "graphics/pointcloud.h"
#include "exporter.h"

class Shader;

struct Fluid {

    double mass;
    double viscosity;
    double gasConstant;
    float colorI; // interface tension color
    float colorS; // surface tension color
    int numParticles;
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

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    // New method to update simulation parameters at runtime
    void updateParameters(
        float new_fluid1_density,
        float new_fluid2_density,
        float new_fluid1_viscosity,
        float new_fluid2_viscosity,
        const Eigen::Vector3d& new_gravity,
        float new_h,
        float new_fluid1_idealGasConstant,
        float new_fluid2_idealGasConstant,
        float new_surfaceTensionThreshold,
        float new_surfaceTensionCoeff,
        float new_interfaceTensionThreshold,
        float new_interfaceTensionCoeff,
        float new_diffusionCoeff
    );

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

private:
    Shape m_shape;
    Shape m_box;
    Shape m_collider;
    QSettings& _settings;
    std::vector<Eigen::Vector3d> m_vertices;
    std::vector<Eigen::Vector3i> m_faces;

    std::vector<Particle *> m_particles;
    std::vector<Fluid*> m_fluids;
    PointCloud m_pointcloud1;
    PointCloud m_pointcloud2;

    Shape m_ground;
    Exporter m_exporter;
    void initGround();
    void initBox();

    // Changed from const to mutable parameters
    float fluid1_density;
    float fluid2_density;
    float fluid1_viscosity;
    float fluid2_viscosity;
    float fluid1_mass;
    float fluid2_mass;
    float fluid1_idealGasConstant;
    float fluid2_idealGasConstant;
    float h;
    float Wpoly6Coeff;
    float Wpoly6GradCoeff;
    float Wpoly6LaplacianCoeff;
    float WspikyGradCoeff;
    float WviscosityLaplacianCoeff;
    Eigen::Vector3d gravity;
    float surfaceTensionThreshold;
    float surfaceTensionCoeff;
    float interfaceTensionThreshold;
    float interfaceTensionCoeff;
    float diffusionCoeff;
    
    // Helper method to update kernel coefficients when h changes
    void updateKernelCoefficients();
};
