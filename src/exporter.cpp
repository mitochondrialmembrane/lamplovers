// ParticleExporter.cpp
#include "Exporter.h"
#include <ImathVec.h>

using namespace Alembic::AbcGeom;

void Exporter::init(const std::string& filename, double fps) {
    m_archive = OArchive(Alembic::AbcCoreOgawa::WriteArchive(), filename);

    OObject topObj = m_archive.getTop();

    m_timeSampling = TimeSamplingPtr(new TimeSampling(1.0 / fps, 0.0));
    m_timeSamplingIndex = m_archive.addTimeSampling(*m_timeSampling);

    m_pointsObj = OPoints(topObj, "particles", m_timeSamplingIndex);
    m_schema = &m_pointsObj.getSchema();
}

void Exporter::addFrame(const std::vector<Eigen::Vector3d>& particles) {
    std::vector<Imath::V3f> positions;
    std::vector<Alembic::Util::uint64_t> ids;
    positions.reserve(particles.size());
    ids.reserve(particles.size());

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        positions.emplace_back(
            static_cast<float>(p.x()),
            static_cast<float>(p.y()),
            static_cast<float>(p.z())
            );
        ids.push_back(static_cast<Alembic::Util::uint64_t>(i));
    }

    Abc::P3fArraySample positionSample(positions);
    Abc::UInt64ArraySample idSample(ids);
    OPointsSchema::Sample sample(positionSample, idSample);

    m_schema->set(sample);
}
