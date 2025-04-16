#include "simulation.h"
#include "graphics/meshloader.h"
#include "util/simulationutils.h"

#include <iostream>
#include <unordered_set>

using namespace Eigen;
using namespace SimulationUtils;

Simulation::Simulation(QSettings& settings) :
    _settings(settings),
    fluid1_density(settings.value("Parameters/fluid1_density").toFloat()),
    fluid1_viscosity(settings.value("Parameters/fluid1_viscosity").toFloat()),
    fluid1_mass(1),
    idealGasConstant(settings.value("Parameters/idealGasConstant").toFloat()),
    surfaceTensionThreshold(settings.value("Parameters/surfaceTensionThreshold").toFloat()),
    surfaceTensionCoeff(settings.value("Parameters/surfaceTensionCoeff").toFloat()),
    gravity(Vector3d(settings.value("Parameters/gravity_x").toFloat(),
                    settings.value("Parameters/gravity_y").toFloat(),
                    settings.value("Parameters/gravity_z").toFloat())),
    h(settings.value("Parameters/smoothingLength").toFloat()),
    m_pointcloud(20)
{
    // Initialize kernel coefficients
    updateKernelCoefficients();
}

// Simulation::Simulation(QSettings& settings) :
//     _settings(settings),
//     fluid1_density(settings.value("Parameters/fluid1_density").toFloat()),
//     fluid1_viscosity(settings.value("Parameters/fluid1_viscosity").toFloat()),
//     fluid1_mass(1),
//     idealGasConstant(settings.value("Parameters/idealGasConstant").toFloat()),
//     surfaceTensionThreshold(settings.value("Parameters/surfaceTensionThreshold").toFloat()),
//     surfaceTensionCoeff(settings.value("Parameters/surfaceTensionCoeff").toFloat()),
//     gravity(Vector3d(settings.value("Parameters/gravity_x").toFloat(),
//                     settings.value("Parameters/gravity_y").toFloat(),
//                     settings.value("Parameters/gravity_z").toFloat())),
//     h(settings.value("Parameters/smoothingLength").toFloat()),
//
//     // Kernel coeffs
//
//     Wpoly6Coeff(315.f / (M_PI * pow(h, 9) * 64)),
//     Wpoly6GradCoeff(945.f / (M_PI * pow(h, 9) * 32)),
//     Wpoly6LaplacianCoeff(945.f / (M_PI * pow(h, 9) * 8)),
//     WspikyGradCoeff(45.f / (M_PI * pow(h, 6))),
//     WviscosityLaplacianCoeff(45.f / (M_PI * pow(h, 6))),    // coeff of spiky kernel grad and viscosity kernel laplacian
//
//     m_pointcloud(20)
//
// {}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    /**if (MeshLoader::loadFaceMesh(":/" + _settings.value("IO/infile").toString().toStdString(), m_vertices, m_faces)) {
        // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
        //    of tetrahedra which are on the exterior surface of the object. Right now, this is
        //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
        //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
        //    are surface faces...

        //m_shape.init(m_vertices, faces, m_tets);
    }**/
    //m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0)));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 20; j++) {
            for (int k = 0; k < 5; k++) {
                m_positions.push_back(Vector3d(0.1 * i + 0.6, 0.1 * j, 0.1 * k - 0.2));
            }
        }
    }

    m_pointcloud.init(m_positions);

    m_particles.resize(m_positions.size());
    for (int i = 0; i < m_positions.size(); i++) {
        m_particles[i] = new Particle{Vector3d(0,0,0), 1, 0, fluid1_density, idealGasConstant, fluid1_viscosity, fluid1_density, 300, 0.5, 1};
    }

    initGround();
    initBox();
}

void Simulation::update(double seconds)
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

    // m_particles[i]->velocity += accelerations[i] * seconds;
    // m_positions[i] += m_particles[i]->velocity * seconds;
    //std::cout << m_positions[i] << std::endl;


    seconds = 0.01;

    std::vector<Vector3d> accelerations;
    accelerations.resize(m_positions.size());
    for (int i = 0; i < m_positions.size(); i++) {
        accelerations[i] = (fPressure(i) + fViscosity(i)) / m_particles[i]->density + gravity;
    }

    std::vector<Vector3d> newAccelerations = accelerations;
    for (int i = 0; i < m_positions.size(); i++) {

        // Leapfrog

        m_positions[i] += m_particles[i]->velocity * seconds + 0.5 * accelerations[i] * seconds * seconds; // Update position
    }

    for (int i = 0; i < m_positions.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(m_particles[i]->density, m_particles[i]->restDensity); // Update density+pressure using new position
    }

    for (int i = 0; i < m_positions.size(); i++) {

        newAccelerations[i] = (fPressure(i) + fViscosity(i) + fSurfaceTension(i)) / m_particles[i]->density + gravity; // Update acceleration
    }

    for (int i = 0; i < m_positions.size(); i++) {

        m_particles[i]->velocity += 0.5 * (accelerations[i] + newAccelerations[i]) * seconds; // Update velocity using both position and acceleration
        evaluateCollisions(i);

    }
    for (int i = 0; i < m_positions.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(m_particles[i]->density, m_particles[i]->restDensity);
    }

    m_pointcloud.setPoints(m_positions);
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
    m_pointcloud.draw(shader);
    m_box.draw(shader);
    //m_collider.draw(shader);
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

void Simulation::initBox()
{
    float a = 1.5;
    float b = 4;
    std::vector<Eigen::Vector3d> vertices = {
        {-a, 0, -a}, {-a, 0, a}, {a, 0, -a}, {a, 0, a},
        {-a, b, -a}, {-a, b, a}, {a, b, -a}, {a, b, a}
    };

    std::vector<Vector3i> faces = {
        // Bottom face (y = 0)
        {0, 2, 3}, {0, 3, 1},
        // Top face (y = b)
        {4, 5, 7}, {4, 7, 6},
        // Front face (z = -a)
        {0, 4, 6}, {0, 6, 2},
        // Back face (z = a)
        {1, 3, 7}, {1, 7, 5},
        // Left face (x = -a)
        {0, 1, 5}, {0, 5, 4},
        // Right face (x = a)
        {2, 6, 7}, {2, 7, 3}
    };
    m_shape.init(vertices, faces);
}

double Simulation::density_S(int i) {
    double out = 0;
    for (int j = 0; j < m_positions.size(); j++) {
        // if (i == j) continue; not necessary??
        double r_squared = (m_positions[i] - m_positions[j]).squaredNorm() + 0.00001;
        if (r_squared < h * h) {
            out += m_particles[j]->mass * Wpoly6Coeff * pow((h * h - r_squared), 3);
        }
    }
    return out;
}

double Simulation::calculatePressure(double density, double restDensity) {
    return idealGasConstant * (density - restDensity);
}

Vector3d Simulation::fPressure(int i) {
    double pressureI = m_particles[i]->pressure;
    Vector3d out(0,0,0);
    for (int j = 0; j < m_positions.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_positions[i] - m_positions[j];
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-WspikyGradCoeff * pow((h - r_norm), 2) / r_norm) * r;
            out += particleJ->mass * (particleJ->pressure + pressureI) / (2 * particleJ->density) * W_Grad;
        }
    }
    return -out;
}

Vector3d Simulation::fViscosity(int i) {
    Vector3d velocityI = m_particles[i]->velocity;
    Vector3d out(0,0,0);
    for (int j = 0; j < m_positions.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_positions[i] - m_positions[j];
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            double W_Laplacian = WviscosityLaplacianCoeff * (h - r_norm);
            out += (particleJ->mass / particleJ->density * W_Laplacian) * (particleJ->velocity - velocityI);
        }
    }

    return m_particles[i]->viscosity * out;
}

Vector3d Simulation::fSurfaceTension(int i) {
    Vector3d normal(0,0,0);
    for (int j = 0; j < m_positions.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_positions[i] - m_positions[j];
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-Wpoly6GradCoeff * pow((h - r_norm), 2)) * r;
            normal += (particleJ->mass / particleJ->density) * W_Grad;
        }
    }
    // std::cout << normal.norm() << std::endl;
    if (normal.norm() < surfaceTensionThreshold) return Vector3d(0,0,0);

    double C_S_Laplacian = 0;
    for (int j = 0; j < m_positions.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_positions[i] - m_positions[j];
        double r_squaredNorm = r.squaredNorm() + 0.00001;
        double h_squared = h * h;
        if (r_squaredNorm < h_squared) {
            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));
            C_S_Laplacian += (particleJ->mass / particleJ->density) * W_Laplacian;
        }
    }
    return (-surfaceTensionCoeff * C_S_Laplacian) * normal.normalized();
}

void Simulation::evaluateCollisions(int i) {
    Vector3d displacement;
    if ((displacement = checkCollision(m_positions[i])) != Vector3d(0,0,0)) {
        m_positions[i] += displacement;
        displacement.normalize();
        m_particles[i]->velocity -= 1.7 * (m_particles[i]->velocity.dot(displacement) * displacement);
    }
}

Vector3d Simulation::checkCollision(Vector3d pos) {
    float wallX = 1.5;
    float wallZ = 0.2;
    float ceiling = 2;

    Vector3d returnVector(0,0,0);
    if (pos[1] < 0) returnVector += Vector3d(0,1,0) * -pos[1];
    if (pos[1] > ceiling) returnVector += Vector3d(0,-1,0) * (pos[1] - ceiling);
    if (pos[0] > wallX) returnVector += Vector3d(-1,0,0) * (pos[0] - wallX);
    if (pos[0] < -wallX) returnVector += Vector3d(1,0,0) * (-wallX - pos[0]);
    if (pos[2] > wallZ) returnVector += Vector3d(0,0,-1) * (pos[2] - wallZ);
    if (pos[2] < -wallZ) returnVector += Vector3d(0,0,1) * (-wallZ - pos[2]);
    return returnVector;
}

void Simulation::updateParameters(
    float new_fluid1_density,
    float new_fluid1_viscosity,
    const Eigen::Vector3d& new_gravity,
    float new_h,
    float new_idealGasConstant,
    float new_surfaceTensionThreshold,
    float new_surfaceTensionCoeff
) {
    // Update parameters
    fluid1_density = new_fluid1_density;
    fluid1_viscosity = new_fluid1_viscosity;
    gravity = new_gravity;
    idealGasConstant = new_idealGasConstant;
    surfaceTensionThreshold = new_surfaceTensionThreshold;
    surfaceTensionCoeff = new_surfaceTensionCoeff;
    
    // Only update h and kernel coefficients if h has changed
    if (h != new_h) {
        h = new_h;
        updateKernelCoefficients();
    }

    // Update rest density, viscosity

    for (int i = 0; i < m_positions.size(); i++) {
        m_particles[i]->restDensity = fluid1_density;
        m_particles[i]->viscosity = fluid1_viscosity;
    }
    
    // Update the density and pressure of all particles
    for (int i = 0; i < m_positions.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(m_particles[i]->density, m_particles[i]->restDensity);
    }
}

void Simulation::updateKernelCoefficients() {
    // Recalculate kernel coefficients based on new h value
    Wpoly6Coeff = 315.f / (M_PI * pow(h, 9) * 64);
    Wpoly6GradCoeff = 945.f / (M_PI * pow(h, 9) * 32);
    Wpoly6LaplacianCoeff = 945.f / (M_PI * pow(h, 9) * 8);
    WspikyGradCoeff = 45.f / (M_PI * pow(h, 6));
    WviscosityLaplacianCoeff = 45.f / (M_PI * pow(h, 6));
}

void Simulation::reinitialize()
{
    // Clear existing particles
    for (auto particle : m_particles) {
        delete particle;
    }
    m_particles.clear();
    m_positions.clear();
    
    // Reinitialize simulation
    init();
}
