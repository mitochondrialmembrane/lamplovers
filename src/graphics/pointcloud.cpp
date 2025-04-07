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

void PointCloud::init(const std::vector<Eigen::Vector3d> &points)
{
    glGenBuffers(1, &m_vbo);
    glGenVertexArrays(1, &m_vao);

    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(double) * points.size() * 3, static_cast<const void *>(points.data()), GL_STATIC_DRAW);

    // Define vertex attributes
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    m_points = points;
    m_numPoints = points.size();

    // Random color
    m_red = static_cast<float>(rand()) / RAND_MAX;
    m_blue = static_cast<float>(rand()) / RAND_MAX;
    m_green = static_cast<float>(rand()) / RAND_MAX;
    m_alpha = 1;
}

void PointCloud::setPoints(const std::vector<Eigen::Vector3d> &points)
{
    if(points.size() != m_numPoints) {
        std::cerr << "You can't set vertices to a vector that is a different length that what point cloud was inited with" << std::endl;
        return;
    }
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(double) * points.size() * 3, points.data());
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
