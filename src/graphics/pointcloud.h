#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <GL/glew.h>
#include <vector>

#include <Eigen/Dense>

class Shader;

struct Point {

    Eigen::Vector3d position;
    double temperature;

};

class PointCloud
{
public:
    PointCloud(float point_size);

    void init(const std::vector<Point> &particles, int i);
    void setPoints(const std::vector<Point> &points);

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

    std::vector<Eigen::Vector3d> m_positions;
    std::vector<double> m_temperatures;

    Eigen::Matrix4f m_modelMatrix;
};

#endif
