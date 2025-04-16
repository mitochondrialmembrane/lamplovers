// ParticleExporter.h
#pragma once

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>
#include <Eigen/Core>
#include <vector>
#include <string>

class Exporter {
public:
    Exporter() = default;

    // Initialize the Alembic archive
    void init(const std::string& filename, double fps = 24.0);

    // Add a frame with particles for both fluids
    void addFrame(const std::vector<Eigen::Vector3d>& fluid1Particles,
                  const std::vector<Eigen::Vector3d>& fluid2Particles);

private:
    Alembic::AbcGeom::OArchive m_archive;    // Alembic archive to store exported data
    Alembic::AbcGeom::OPoints m_fluid1PointsObj;  // Points object for fluid 1
    Alembic::AbcGeom::OPoints m_fluid2PointsObj;  // Points object for fluid 2
    Alembic::AbcGeom::OPointsSchema* m_fluid1Schema = nullptr;  // Points schema for fluid 1
    Alembic::AbcGeom::OPointsSchema* m_fluid2Schema = nullptr;  // Points schema for fluid 2
    Alembic::AbcGeom::TimeSamplingPtr m_timeSampling;  // Time sampling information
    uint32_t m_timeSamplingIndex = 0;  // Time sampling index
};
