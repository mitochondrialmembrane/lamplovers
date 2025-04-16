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

    // Add a frame at the given time
    void addFrame(const std::vector<Eigen::Vector3d>& particles);

private:
    Alembic::AbcGeom::OArchive m_archive;
    Alembic::AbcGeom::OPoints m_pointsObj;
    Alembic::AbcGeom::OPointsSchema* m_schema = nullptr;
    Alembic::AbcGeom::TimeSamplingPtr m_timeSampling;
    uint32_t m_timeSamplingIndex = 0;
};
