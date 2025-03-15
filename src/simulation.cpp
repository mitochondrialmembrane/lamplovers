#include "simulation.h"
#include "graphics/meshloader.h"
#include "util/simulationutils.h"

#include <iostream>
#include <unordered_set>

using namespace Eigen;
using namespace SimulationUtils;

Simulation::Simulation(QSettings& settings) :
    _settings(settings),
    gravity(Vector3d(settings.value("Parameters/gravityx").toFloat(),
                    settings.value("Parameters/gravityy").toFloat(),
                    settings.value("Parameters/gravityx").toFloat())),
    lambda(settings.value("Parameters/lambda").toFloat()),
    mu(settings.value("Parameters/mu").toFloat()),
    phi(settings.value("Parameters/phi").toFloat()),
    psi(settings.value("Parameters/psi").toFloat()),
    rho(settings.value("Parameters/rho").toFloat()),
    cor(settings.value("Parameters/cor").toFloat()),
    friction(settings.value("Parameters/friction").toFloat())
{}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    std::unordered_map<Vector3i, Vector3i, FaceHash> faceSet;
    if (MeshLoader::loadTetMesh(":/" + _settings.value("IO/infile").toString().toStdString(), m_vertices, m_tets)) {
        // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
        //    of tetrahedra which are on the exterior surface of the object. Right now, this is
        //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
        //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
        //    are surface faces...
        std::vector<Vector3i> faces;
        m_material_coords = m_vertices;
        v_masses.resize(m_vertices.size());
        v_velocities.resize(m_vertices.size());
        for (int i = 0; i < m_vertices.size(); i++) {
            v_velocities[i] = Vector3d(0,0,0);
            v_masses[i] = 0;
        }

        Vector4<Vector3i> indices(Vector3i(1,2,3), Vector3i(0,3,2), Vector3i(0,1,3), Vector3i(0,2,1));
        m_betas.resize(m_tets.size());
        m_volumes.resize(m_tets.size());
        t_normals.resize(m_tets.size());
        t_areas.resize(m_tets.size());

        for (int i = 0; i < m_tets.size(); i++) {
            Vector4i tet = m_tets[i];
            m_volumes[i] = tetrahedronVolume(m_vertices[tet[0]], m_vertices[tet[1]], m_vertices[tet[2]], m_vertices[tet[3]]);
            double m = m_volumes[i] * rho / 4.f;

            Matrix4d mat;
            Vector4<Vector3d> normals;
            Vector4<double> areas;

            for (int j = 0; j < 4; j++) {
                v_masses[tet[j]] += m;

                // add to face surface set
                Vector3i index = indices[j];
                std::array<int, 3> values = {tet[index[0]], tet[index[1]], tet[index[2]]};
                Vector3i face(values[0], values[1], values[2]);
                std::sort(values.begin(), values.end());
                Vector3i faceSorted(values[0], values[1], values[2]);

                if (faceSet.contains(faceSorted)) {
                    faceSet.erase(faceSorted);
                }
                else {
                    faceSet[faceSorted] = face;
                }

                // initialize beta matrix
                Vector3d vertex(m_vertices[m_tets[i][j]]);
                mat.col(j) = Vector4d(vertex[0], vertex[1], vertex[2], 1);
                Vector3d normal = (m_vertices[tet[index[1]]] - m_vertices[tet[index[0]]]).cross(m_vertices[tet[index[2]]] - m_vertices[tet[index[0]]]);
                normals[j] = normal.normalized();
                areas[j] = normal.norm() * 0.5;
            }
            m_betas[i] = mat.inverse();
            t_normals[i] = normals;
            t_areas[i] = areas;
        }

        for (const auto& [key, value] : faceSet) {
            faces.push_back(value);
        }

        m_shape.init(m_vertices, faces, m_tets);
    }
    m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0)));

    initGround();
    initCollider();
}

void Simulation::update(double seconds, bool paused, Vector2d dragChange, Vector3f look)
{
    // STUDENTS: This method should contain all the time-stepping logic for your simulation.
    //   Specifically, the code you write here should compute new, updated vertex positions for your
    //   simulation mesh, and it should then call m_shape.setVertices to update the display with those
    //   newly-updated vertices.

    // STUDENTS: As currently written, the program will just continually compute simulation timesteps as long
    //    as the program is running (see View::tick in view.cpp) . You might want to e.g. add a hotkey for pausing
    //    the simulation, and perhaps start the simulation out in a paused state.

    // Note that the "seconds" parameter represents the amount of time that has passed since
    // the last update
    if (paused) {
        return;
    }

    std::vector<Vector3d> internalForces;
    std::vector<Vector3d> midpoints;
    internalForces.resize(m_vertices.size());
    midpoints.resize(m_vertices.size());
    std::fill(internalForces.begin(), internalForces.end(), Vector3d(0,0,0));
    Vector4<Vector3i> indices(Vector3i(1,2,3), Vector3i(0,3,2), Vector3i(0,1,3), Vector3i(0,2,1));

    seconds = 0.001;
    for (int i = 0; i < m_vertices.size(); i++) {
        midpoints[i] = m_vertices[i] + v_velocities[i] * seconds / 2;
    }

    for (int i = 0; i < m_tets.size(); i++) {
        Vector4i tet = m_tets[i];
        Matrix4d beta = m_betas[i];

        Matrix<double, 3, 4> P;
        Matrix<double, 3, 4> V;
        // calculate P and V matrices
        for (int j = 0; j < 4; j++) {
            Vector3d point = midpoints[tet[j]];
            Vector3d velocity = v_velocities[tet[j]];

            for (int k = 0; k < 3; k++) {
                P(k, j) = point[k];
                V(k, j) = velocity[k];
            }
        }

        // calculate partials
        Matrix<double, 3, 4> xPartials = P * beta;
        Matrix<double, 3, 4> xDotPartials = V * beta;

        Matrix3d F;
        for (int i = 0; i < 3; i++) {
            F.col(i) = xPartials.col(i);
        }

        Matrix3d Fdot;
        for (int i = 0; i < 3; i++) {
            Fdot.col(i) = xDotPartials.col(i);
        }

        // calculate tensors
        Matrix3d greensStrain = F.transpose() * F - Matrix3d::Identity();
        Matrix3d strainRate = Fdot.transpose() * F + F.transpose() * Fdot;


        Matrix3d elasticStress = lambda * greensStrain.trace() * Matrix3d::Identity() + 2 * mu * greensStrain;
        Matrix3d viscousStress = phi * strainRate.trace() * Matrix3d::Identity() + 2 * psi * strainRate;

        Matrix3d totalStress = elasticStress + viscousStress;

        Vector4<Vector3d> normals = t_normals[i];
        Vector4<double> areas = t_areas[i];

        for (int j = 0; j < 4; j++) {
            Vector3d force = areas[j] * F * (totalStress * normals[j]);
            Vector3i index = indices[j];

            for (int k = 0; k < 3; k++) {
                internalForces[tet[index[k]]] -= force;
            }
        }

    }


    for (int i = 0; i < m_vertices.size(); i++) {
        Vector3d midpointPos = midpoints[i];

        Vector3d force = calculateExternalForces(i, midpointPos, dragChange, look) + internalForces[i];

        v_velocities[i] += force / v_masses[i] * seconds;

        m_vertices[i] += v_velocities[i] * seconds;

        Vector3d collisionNormal;
        if ((collisionNormal = checkCollision(m_vertices[i])).norm() > 1e-6) {
            m_vertices[i] += collisionNormal;
            collisionNormal.normalize();
            Vector3d v_n = -(v_velocities[i].dot(collisionNormal) / collisionNormal.dot(collisionNormal)) * collisionNormal;
            Vector3d v_t = v_velocities[i] + v_n;
            v_n.normalize();
            v_t.normalize();

            v_velocities[i] = v_n * cor + v_t * friction;
        }
    }
    m_shape.setVertices(m_vertices);
}

Vector3d Simulation::calculateExternalForces(int i, Vector3d pos, Vector2d dragChange, Vector3f look) {
    Vector3d force = gravity * v_masses[i];

    Vector3d lookDouble(look[0], look[1], look[2]);

    if (dragChange != Vector2d(0,0)) {
        Vector3d x = lookDouble.cross(Vector3d(0,1,0)).normalized();
        Vector3d y = x.cross(lookDouble).normalized();
        force += (x * dragChange[0] + y * dragChange[1]) * v_masses[i] * 2;
    }

    return force;
}

Vector3d Simulation::checkCollision(Eigen::Vector3d pos) {
    pos[1] += 2;
    if (pos[1] <= 0) {
        return Vector3d(0,-pos[1],0);
    }
    double colliderRad = 1.f;
    double r_2 = pow(pos[0] - 0.5, 2) + pow(pos[1], 2) + pow(pos[2], 2);
    if (r_2 < colliderRad * colliderRad) {
        return -(pos - Vector3d(0.5, 2, 0)).normalized() * (colliderRad - sqrt(r_2)) * 0.5;
    }
    return Vector3d(0,0,0);
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
    m_collider.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

void Simulation::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces);
}

void Simulation::initCollider()
{
    std::unordered_map<Vector3i, Vector3i, FaceHash> faceSet;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector4i> tets;

    if (MeshLoader::loadTetMesh(":/example-meshes/sphere.mesh", vertices, tets)) {
        std::vector<Vector3i> faces;

        Vector4<Vector3i> indices(Vector3i(1,2,3), Vector3i(0,3,2), Vector3i(0,1,3), Vector3i(0,2,1));

        for (int i = 0; i < tets.size(); i++) {
            Vector4i tet = tets[i];

            for (int j = 0; j < 4; j++) {

                // add to face surface set
                Vector3i index = indices[j];
                std::array<int, 3> values = {tet[index[0]], tet[index[1]], tet[index[2]]};
                Vector3i face(values[0], values[1], values[2]);
                std::sort(values.begin(), values.end());
                Vector3i faceSorted(values[0], values[1], values[2]);

                if (faceSet.contains(faceSorted)) {
                    faceSet.erase(faceSorted);
                }
                else {
                    faceSet[faceSorted] = face;
                }
            }
        }

        for (const auto& [key, value] : faceSet) {
            faces.push_back(value);
        }

        m_collider.init(vertices, faces, tets);
    }
    m_collider.setModelMatrix(Affine3f(Eigen::Translation3f(0.5, 0, 0)));
}
