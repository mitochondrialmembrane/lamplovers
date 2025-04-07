#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <GL/glew.h>
#include <vector>

#include <Eigen/Dense>

class Shader;

class PointCloud
{
public:
    PointCloud(float point_size);

    void init(const std::vector<Eigen::Vector3d> &points);
    void setPoints(const std::vector<Eigen::Vector3d> &points);

    void setModelMatrix(const Eigen::Affine3f &model);

    void draw(Shader *shader);

private:
    GLuint m_vao;
    GLuint m_vbo;

    unsigned int m_numPoints;
    float m_pointSize;

    float m_red;
    float m_blue;
    float m_green;
    float m_alpha;

    std::vector<Eigen::Vector3d> m_points;

    Eigen::Matrix4f m_modelMatrix;
};

#endif
