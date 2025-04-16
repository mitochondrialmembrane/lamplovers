#include "Exporter.h"
#include <ImathVec.h>

using namespace Alembic::AbcGeom;

void Exporter::init(const std::string& filename, double fps) {
    m_archive = OArchive(Alembic::AbcCoreOgawa::WriteArchive(), filename);

    OObject topObj = m_archive.getTop();

    m_timeSampling = TimeSamplingPtr(new TimeSampling(1.0 / fps, 0.0));
    m_timeSamplingIndex = m_archive.addTimeSampling(*m_timeSampling);

    // Create two separate OPoints for two fluids
    m_fluid1PointsObj = OPoints(topObj, "fluid1_particles", m_timeSamplingIndex);
    m_fluid2PointsObj = OPoints(topObj, "fluid2_particles", m_timeSamplingIndex);

    // Create schemas for both fluids
    m_fluid1Schema = &m_fluid1PointsObj.getSchema();
    m_fluid2Schema = &m_fluid2PointsObj.getSchema();
}

void Exporter::addFrame(const std::vector<Eigen::Vector3d>& particlesFluid1, const std::vector<Eigen::Vector3d>& particlesFluid2) {
    // Handle Fluid 1 (First Fluid)
    std::vector<Imath::V3f> positionsFluid1;
    std::vector<Alembic::Util::uint64_t> idsFluid1;
    positionsFluid1.reserve(particlesFluid1.size());
    idsFluid1.reserve(particlesFluid1.size());

    for (size_t i = 0; i < particlesFluid1.size(); ++i) {
        const auto& p = particlesFluid1[i];
        positionsFluid1.emplace_back(
            static_cast<float>(p.x()),
            static_cast<float>(p.y()),
            static_cast<float>(p.z())
            );
        idsFluid1.push_back(static_cast<Alembic::Util::uint64_t>(i));
    }

    // Create a sample for Fluid 1
    Abc::P3fArraySample positionSampleFluid1(positionsFluid1);
    Abc::UInt64ArraySample idSampleFluid1(idsFluid1);
    OPointsSchema::Sample sampleFluid1(positionSampleFluid1, idSampleFluid1);

    // Set the sample for Fluid 1
    m_fluid1Schema->set(sampleFluid1);

    // Handle Fluid 2 (Second Fluid)
    std::vector<Imath::V3f> positionsFluid2;
    std::vector<Alembic::Util::uint64_t> idsFluid2;
    positionsFluid2.reserve(particlesFluid2.size());
    idsFluid2.reserve(particlesFluid2.size());

    for (size_t i = 0; i < particlesFluid2.size(); ++i) {
        const auto& p = particlesFluid2[i];
        positionsFluid2.emplace_back(
            static_cast<float>(p.x()),
            static_cast<float>(p.y()),
            static_cast<float>(p.z())
            );
        idsFluid2.push_back(static_cast<Alembic::Util::uint64_t>(i));
    }

    // Create a sample for Fluid 2
    Abc::P3fArraySample positionSampleFluid2(positionsFluid2);
    Abc::UInt64ArraySample idSampleFluid2(idsFluid2);
    OPointsSchema::Sample sampleFluid2(positionSampleFluid2, idSampleFluid2);

    // Set the sample for Fluid 2
    m_fluid2Schema->set(sampleFluid2);
}
