#pragma once

#include "graphics/shape.h"
#include <QtCore>

class Shader;

class Simulation
{
public:
    Simulation(QSettings& settings);

    void init();

    void update(double seconds, bool paused, Eigen::Vector2d dragChange, Eigen::Vector3f look);

    void draw(Shader *shader);

    void toggleWire();

    struct FaceHash {
        std::size_t operator()(Eigen::Vector3i vertices) const {
            // Hash the sorted vertex indices
            return std::hash<int>()(vertices[0]) ^
                   (std::hash<int>()(vertices[1]) << 1) ^
                   (std::hash<int>()(vertices[2]) << 2);
        }
    };
private:
    Shape m_shape;
    Shape m_collider;
    QSettings& _settings;
    std::vector<Eigen::Vector3d> v_velocities;
    std::vector<double> v_masses;
    std::vector<Eigen::Vector3d> m_vertices;
    std::vector<Eigen::Vector4i> m_tets;
    std::vector<double> m_volumes;
    std::vector<Eigen::Vector3d> m_material_coords;
    std::vector<Eigen::Vector4<Eigen::Vector3d>> t_normals;
    std::vector<Eigen::Vector4<double>> t_areas;
    std::vector<Eigen::Matrix4d> m_betas;

    Shape m_ground;
    void initGround();
    void initCollider();
    Eigen::Vector3d calculateExternalForces(int i, Eigen::Vector3d pos, Eigen::Vector2d dragChange, Eigen::Vector3f look);
    Eigen::Vector3d checkCollision(Eigen::Vector3d pos);

    Eigen::Vector3d gravity;
    const float lambda;
    const float mu;
    const float phi;
    const float psi;
    const float rho;
    const float cor;
    const float friction;
};
