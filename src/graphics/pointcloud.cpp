#include "pointcloud.h"

#include <iostream>

#include "graphics/shader.h"

using namespace Eigen;

PointCloud::PointCloud(float point_size)
    : m_vao(-1),
    m_numPoints(),
    m_pointSize(point_size),
    m_modelMatrix(Eigen::Matrix4f::Identity())
{
}

void PointCloud::init(const std::vector<Point>& points, int i)
{

    std::vector<Vector3d> positions;
    std::vector<double> temperatures;

    for (int i = 0; i < points.size(); i++) {

        positions.push_back(points[i].position);
        temperatures.push_back(points[i].temperature);

    }

    glGenBuffers(1, &m_vbo);
    glGenVertexArrays(1, &m_vao);

    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(double) * positions.size() * sizeof(Point), static_cast<const void *>(points.data()), GL_STATIC_DRAW);

    // Define vertex attributes
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(Point), 0);
    glEnableVertexAttribArray(0);


    glVertexAttribPointer(
        1,
        1, GL_DOUBLE,
        GL_FALSE,
        sizeof(Point),
        (void*)offsetof(Point, temperature)
        );
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    m_positions = positions;
    m_temperatures = temperatures;
    m_numPoints = points.size();

    // Random color
    m_red = 0;
    m_blue = 1 - i;
    m_green = i;
    m_alpha = i * 0.75 + 0.25;
}

void PointCloud::setPoints(const std::vector<Point> &points)
{
    if(points.size() != m_numPoints) {
        std::cerr << "You can't set vertices to a vector that is a different length that what point cloud was inited with" << std::endl;
        return;
    }

    std::vector<Vector3d> positions;
    std::vector<double> temperatures;

    for (int i = 0; i < points.size(); i++) {

        positions.push_back(points[i].position);
        temperatures.push_back(points[i].temperature);

    }

    m_positions = positions;
    m_temperatures = temperatures;

    // for (int i = 0; i < m_temperatures.size(); i++) {

    //     std::cout << m_temperatures[i] << std::endl;

    // }

    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(double) * points.size() * 4, points.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void PointCloud::setModelMatrix(const Eigen::Affine3f &model)
{
    m_modelMatrix = model.matrix();
}

void PointCloud::draw(Shader *shader)
{
    shader->setUniform("model", m_modelMatrix);
    shader->setUniform("red",   m_red);
    shader->setUniform("green", m_green);
    shader->setUniform("blue",  m_blue);
    shader->setUniform("alpha", m_alpha);
    shader->setUniform("pointSize", m_pointSize);
    shader->setUniform("u_renderMode", 1);

    glBindVertexArray(m_vao);
    glDrawArrays(GL_POINTS, 0, m_numPoints);  // ✅ Use glDrawArrays for point clouds
    glBindVertexArray(0);
}
