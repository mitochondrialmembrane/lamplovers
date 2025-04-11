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
};

class Simulation
{
public:
    Simulation(QSettings& settings);

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    double density_S(int i);
    double calculatePressure(double density);
    Eigen::Vector3d fPressure(int i);
    Eigen::Vector3d fViscosity(int i);
    Eigen::Vector3d fSurfaceTension(int i);
    Eigen::Vector3d checkCollision(Eigen::Vector3d pos);
    void evaluateCollisions(int i);

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

    const float fluid1_density;
    const float fluid1_viscosity;
    const float fluid1_mass;
    const float idealGasConstant;
    const float h;
    const float Wpoly6Coeff;
    const float Wpoly6GradCoeff;
    const float Wpoly6LaplacianCoeff;
    const float WspikyGradCoeff;
    const float WviscosityLaplacianCoeff;
    const Eigen::Vector3d gravity;
    const double surfaceTensionThreshold;
    const double surfaceTensionCoeff;
};
