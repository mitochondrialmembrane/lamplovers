#pragma once

#include "graphics/shape.h"
#include <QtCore>
#include "graphics/pointcloud.h"

class Shader;

struct Particle {
    Eigen::Vector3d velocity;
    double mass;
    double pressure;
    double density;

    // NEW ATTRIBUTES FOR MULTIPLE FLUIDS
    double gasConstant;
    double viscosity;
    double restDensity;
    double temperature; // in celsius
    float colorI; // interface tension color
    float colorS; // surface tension color

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
        float fluid1_density,
        float fluid1_viscosity,
        const Eigen::Vector3d& gravity,
        float h,
        float idealGasConstant,
        float surfaceTensionThreshold,
        float surfaceTensionCoeff
    );

    double density_S(int i);
    double calculatePressure(double density, double restDensity);
    Eigen::Vector3d fPressure(int i);
    Eigen::Vector3d fViscosity(int i);
    Eigen::Vector3d fSurfaceTension(int i);
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

    std::vector<Eigen::Vector3d> m_positions;
    std::vector<Particle *> m_particles;
    PointCloud m_pointcloud;

    Shape m_ground;
    void initGround();
    void initBox();

    // Changed from const to mutable parameters
    float fluid1_density;
    float fluid1_viscosity;
    float fluid1_mass;
    float idealGasConstant;
    float h;
    float Wpoly6Coeff;
    float Wpoly6GradCoeff;
    float Wpoly6LaplacianCoeff;
    float WspikyGradCoeff;
    float WviscosityLaplacianCoeff;
    Eigen::Vector3d gravity;
    double surfaceTensionThreshold;
    double surfaceTensionCoeff;
    
    // Helper method to update kernel coefficients when h changes
    void updateKernelCoefficients();
};
