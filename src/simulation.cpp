#include "simulation.h"
#include "graphics/meshloader.h"
#include "util/simulationutils.h"

#include <iostream>
#include <unordered_set>
#include <omp.h>

using namespace Eigen;
using namespace SimulationUtils;

Simulation::Simulation(QSettings& settings) :
    _settings(settings),
    fluid1_density(settings.value("Parameters/fluid1_density").toFloat()),
    fluid2_density(settings.value("Parameters/fluid2_density").toFloat()),
    fluid1_viscosity(settings.value("Parameters/fluid1_viscosity").toFloat()),
    fluid2_viscosity(settings.value("Parameters/fluid2_viscosity").toFloat()),
    fluid1_mass(1),
    fluid2_mass(2),
    fluid1_idealGasConstant(settings.value("Parameters/fluid1_idealGasConstant").toFloat()),
    fluid2_idealGasConstant(settings.value("Parameters/fluid2_idealGasConstant").toFloat()),
    surfaceTensionThreshold(settings.value("Parameters/surfaceTensionThreshold").toFloat()),
    surfaceTensionCoeff(settings.value("Parameters/surfaceTensionCoeff").toFloat()),
    interfaceTensionThreshold(settings.value("Parameters/interfaceTensionThreshold").toFloat()),
    interfaceTensionCoeff(settings.value("Parameters/interfaceTensionCoeff").toFloat()),
    diffusionCoeff(settings.value("Parameters/diffusionCoeff").toFloat()),
    gravity(Vector3d(settings.value("Parameters/gravity_x").toFloat(),
                     settings.value("Parameters/gravity_y").toFloat(),
                     settings.value("Parameters/gravity_z").toFloat())),
    h(settings.value("Parameters/smoothingLength").toFloat()),
    m_pointcloud1(20),
    m_pointcloud2(20),
    radius (0.5),
    ceiling (1.5),
    coneTop (3)
{
    // Initialize kernel coefficients
    updateKernelCoefficients();
    m_exporter.init("test.abc", 100);
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

    // make two fluids

    Fluid* fluid1 = new Fluid{fluid1_mass, fluid1_viscosity, fluid1_idealGasConstant, 0.5, 1, 0};
    Fluid* fluid2 = new Fluid{fluid2_mass, fluid2_viscosity, fluid2_idealGasConstant, -0.5, 1, 0};

    m_fluids.push_back(fluid1);
    m_fluids.push_back(fluid2);

    // initialize two sets of particles

    std::vector<Vector3d> fluid1Positions;
    std::vector<Vector3d> fluid2Positions;
    int d = 15; // Change to increase the density of particles

    #pragma omp parallel
    {
        std::vector<Vector3d> local_fluid1Positions;
        std::vector<Vector3d> local_fluid2Positions;
        
        #pragma omp for collapse(3)
        for (int i = -radius * d; i < radius * d; i++) {
            for (int j = -radius * d; j < radius * d; j++) {
                for (int k = 0; k < ceiling * d; k++) {
                    Vector3d pos((float) i / d, (float) k / d, (float) j / d);
                    if (checkCollision(pos) == Vector3d(0,0,0)) {
                        if (k > 0.3 * d) local_fluid1Positions.push_back(pos);
                        else local_fluid2Positions.push_back(pos);
                    }
                }
            }
        }
        
        // Merge thread-local results to the shared vectors
        #pragma omp critical
        {
            fluid1Positions.insert(fluid1Positions.end(), local_fluid1Positions.begin(), local_fluid1Positions.end());
            fluid2Positions.insert(fluid2Positions.end(), local_fluid2Positions.begin(), local_fluid2Positions.end());
        }
    }
    /**
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 5; k++) {
                fluid2Positions.push_back(Vector3d(0.1 * i - 0.5, 0.1 * j, 0.1 * k - 0.2));
            }
        }
    }**/

    // for (int i = 0; i < 10; i++) {
    //     for (int j = 0; j < 10; j++) {
    //         for (int k = 0; k < 5; k++) {
    //            fluid1Positions.push_back(Vector3d(0.1 * i + 0.6, 0.1 * j, 0.1 * k - 0.2));
    //         }
    //     }
    // }

    // for (int i = 0; i < 10; i++) {
    //     for (int j = 0; j < 10; j++) {
    //         for (int k = 0; k < 5; k++) {
    //             fluid2Positions.push_back(Vector3d(0.1 * i - 1.5, 0.1 * j, 0.1 * k - 0.2));
    //         }
    //     }
    // }

    std::vector<Point> fluid1Points;
    std::vector<Point> fluid2Points;

    for (int i = 0; i < fluid1Positions.size(); i++) {
        Particle* newParticle = new Particle{fluid1Positions[i], Vector3d(0, 0, 0), 0, fluid1_density, fluid1_density, 300, fluid1};

        m_particles.push_back(newParticle);

        fluid1Points.push_back(Point{fluid1Positions[i], 300});
    }

    for (int i = 0; i < fluid2Positions.size(); i++) {
        Particle* newParticle = new Particle{fluid2Positions[i], Vector3d(0, 0, 0), 0, fluid2_density, fluid2_density, 300, fluid2};
        m_particles.push_back(newParticle);

        fluid2Points.push_back(Point{fluid2Positions[i], 300});
    }

    m_fluids[0]->numParticles = fluid1Positions.size();
    m_fluids[1]->numParticles = fluid2Positions.size();



    m_pointcloud1.init(fluid1Points, 0);
    m_pointcloud2.init(fluid2Points, 1);

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

    // Slowly heat the bottom + cool the top

    for (int i = 0; i < m_particles.size(); i++) {

        if (m_particles[i]->position[1] < 0.1) {

            m_particles[i]->temperature += 0.1;

        }

        if (m_particles[i]->position[1] > 0.9) {

            m_particles[i]->temperature -= 0.1;

        }

    }

    // Heat equation for diffusion using Euler

    float TEMPERATURE_DENSITY_CONSTANT = 600000;

    std::vector<double> deltaTValues;
    deltaTValues.resize(m_particles.size());

#pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        double deltaT = calculateTemperatureDiffusionStep(i);
        deltaTValues[i] = deltaT;
    }

#pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->temperature += deltaTValues[i] * seconds;
    }

    // Only liquid 2 (on the bottom) is temperature dependent?

    for (int i = m_fluids[0]->numParticles; i < m_particles.size(); i++) {

        m_particles[i]->restDensity = TEMPERATURE_DENSITY_CONSTANT / m_particles[i]->temperature;

    }

    std::vector<Vector3d> accelerations;
    accelerations.resize(m_particles.size());

#pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        accelerations[i] = (fPressure(i) + fViscosity(i)) / m_particles[i]->density + gravity;
    }

    std::vector<Vector3d> newAccelerations = accelerations;
#pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        // Leapfrog
        m_particles[i]->position += m_particles[i]->velocity * seconds + 0.5 * accelerations[i] * seconds * seconds; // Update position
    }

#pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(i, m_particles[i]->density, m_particles[i]->restDensity); // Update density+pressure using new position
    }

#pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {

        newAccelerations[i] = (fPressure(i) + fViscosity(i) + fSurfaceTension(i) + fInterfaceTension(i)) / m_particles[i]->density + gravity; // Update acceleration
    }

    for (int i = 0; i < m_particles.size(); i++) {

        m_particles[i]->velocity += 0.5 * (accelerations[i] + newAccelerations[i]) * seconds; // Update velocity using both position and acceleration
        evaluateCollisions(i);

    }
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(i, m_particles[i]->density, m_particles[i]->restDensity);
    }

    std::vector<Point> fluid1Points;
    std::vector<Point> fluid2Points;
    std::vector<Vector3d> fluid1Positions;
    std::vector<Vector3d> fluid2Positions;

    for (int i = 0; i < m_fluids[0]->numParticles; i++) {

        fluid1Points.push_back(Point{m_particles[i]->position, m_particles[i]->temperature});
        fluid1Positions.push_back(m_particles[i]->position);

    }

    for (int i = 0; i < m_fluids[1]->numParticles; i++) {

        fluid2Points.push_back(Point{m_particles[m_fluids[0]->numParticles + i]->position, m_particles[m_fluids[0]->numParticles + i]->temperature});
        fluid2Positions.push_back(m_particles[m_fluids[0]->numParticles + i]->position);

    }

    m_pointcloud1.setPoints(fluid1Points);
    m_pointcloud2.setPoints(fluid2Points);

    m_exporter.addFrame(fluid1Positions, fluid2Positions);
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
    m_pointcloud1.draw(shader);
    m_pointcloud2.draw(shader);
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

    // CHANGE FOR MULTI FLUIDS: Still iterate over all particles?? (Check later)

    double out = 0;
    // Using reduction to safely accumulate the sum in parallel
    #pragma omp parallel for reduction(+:out)
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue; not necessary??
        //
        // Each thread calculates contributions from some particles
        // and adds them to its private copy of 'out'
        double r_squared = (m_particles[i]->position - m_particles[j]->position).squaredNorm() + 0.00001;
        if (r_squared < h * h) {
            out += m_particles[j]->fluid->mass * Wpoly6Coeff * pow((h * h - r_squared), 3);
        }
    }

    // At the end, all private copies are added together into the shared 'out'
    return out;
}

double Simulation::calculateTemperatureDiffusionStep(int i) {

    Particle *particleI = m_particles[i];

    double out = 0;

    for (int j = 0; j < m_particles.size(); j++) {

        Particle *particleJ = m_particles[j];
        Vector3d r = particleI->position - particleJ->position;
        double r_squaredNorm = r.squaredNorm() + 0.00001;
        double h_squared = h * h;
        if (r_squaredNorm < h_squared) {

            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));

            out += diffusionCoeff * (particleJ->fluid->mass / particleJ->density * W_Laplacian) * (particleJ->temperature - particleI->temperature);

        }

    }

    return out;

}

double Simulation::calculatePressure(int i, double density, double restDensity) {
    return m_particles[i]->fluid->gasConstant * (density - restDensity);
}

Vector3d Simulation::fPressure(int i) {
    double pressureI = m_particles[i]->pressure;
    Vector3d out(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-WspikyGradCoeff * pow((h - r_norm), 2) / r_norm) * r;
            out += particleJ->fluid->mass * (particleJ->pressure + pressureI) / (2 * particleJ->density) * W_Grad;
        }
    }
    return -out;
}

Vector3d Simulation::fViscosity(int i) {
    Particle *particleI = m_particles[i];
    Vector3d out(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            double W_Laplacian = WviscosityLaplacianCoeff * (h - r_norm);
            double avgViscosity = (particleJ->fluid->viscosity + particleI->fluid->viscosity) / 2.0;
            out += avgViscosity * (particleJ->fluid->mass / particleJ->density * W_Laplacian) * (particleJ->velocity - particleI->velocity);
        }
    }

    return out;
}

Vector3d Simulation::fSurfaceTension(int i) {
    Vector3d normal(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-Wpoly6GradCoeff * pow((h - r_norm), 2)) * r;
            normal += (particleJ->fluid->mass / particleJ->density) * W_Grad;
        }
    }
    // std::cout << normal.norm() << std::endl;
    if (normal.norm() < surfaceTensionThreshold) return Vector3d(0,0,0);

    double C_S_Laplacian = 0;
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_squaredNorm = r.squaredNorm() + 0.00001;
        double h_squared = h * h;
        if (r_squaredNorm < h_squared) {
            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));
            C_S_Laplacian += (particleJ->fluid->mass / particleJ->density) * W_Laplacian;
        }
    }
    return (-surfaceTensionCoeff * C_S_Laplacian) * normal.normalized();
}

Vector3d Simulation::fInterfaceTension(int i) {
    Vector3d normal(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-Wpoly6GradCoeff * pow((h - r_norm), 2)) * r;
            normal += particleJ->fluid->colorI * (particleJ->fluid->mass / particleJ->density) * W_Grad;
        }
    }
    // std::cout << normal.norm() << std::endl;
    if (normal.norm() < interfaceTensionThreshold) return Vector3d(0,0,0);

    double C_I_Laplacian = 0;
    for (int j = 0; j < m_particles.size(); j++) {
        // if (i == j) continue;
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_squaredNorm = r.squaredNorm() + 0.00001;
        double h_squared = h * h;
        if (r_squaredNorm < h_squared) {
            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));
            C_I_Laplacian += particleJ->fluid->colorI * (particleJ->fluid->mass / particleJ->density) * W_Laplacian;
        }
    }
    return (-interfaceTensionCoeff * C_I_Laplacian) * normal.normalized();
}

void Simulation::evaluateCollisions(int i) {
    Vector3d displacement;
    if ((displacement = checkCollision(m_particles[i]->position)) != Vector3d(0,0,0)) {
        m_particles[i]->position += displacement;
        displacement.normalize();
        m_particles[i]->velocity -= 1.7 * (m_particles[i]->velocity.dot(displacement) * displacement);
    }
}

Vector3d Simulation::checkCollision(Vector3d pos) {
    // float wallX = 1.5;
    // float wallZ = 0.2;
    // float ceiling = 2;

    /**float wallX = 0.5;
    float wallZ = 0.2;
    float ceiling = 1;

    Vector3d returnVector(0,0,0);
    if (pos[1] < 0) returnVector += Vector3d(0,1,0) * -pos[1];
    if (pos[1] > ceiling) returnVector += Vector3d(0,-1,0) * (pos[1] - ceiling);
    if (pos[0] > wallX) returnVector += Vector3d(-1,0,0) * (pos[0] - wallX);
    if (pos[0] < -wallX) returnVector += Vector3d(1,0,0) * (-wallX - pos[0]);
    if (pos[2] > wallZ) returnVector += Vector3d(0,0,-1) * (pos[2] - wallZ);
    if (pos[2] < -wallZ) returnVector += Vector3d(0,0,1) * (-wallZ - pos[2]);
    return returnVector;**/

    float r = sqrt(pos[0] * pos[0] + pos[2] * pos[2]);
    float radiusAtY = (coneTop - pos[1]) / coneTop * radius;
    Vector3d returnVector(0,0,0);
    if (pos[1] < 0) returnVector += Vector3d(0,1,0) * -pos[1];
    if (pos[1] > ceiling) returnVector += Vector3d(0,-1,0) * (pos[1] - ceiling);
    if (r > radiusAtY) returnVector += (Vector3d(-pos[0], 0, -pos[2]).normalized() * coneTop - Vector3d(0, -radius, 0)).normalized() * (r - radiusAtY);
    return returnVector;
}



void Simulation::updateParameters(
    float new_fluid1_density,
    float new_fluid2_density,
    float new_fluid1_viscosity,
    float new_fluid2_viscosity,
    const Eigen::Vector3d& new_gravity,
    float new_h,
    float new_fluid1_idealGasConstant,
    float new_fluid2_idealGasConstant,
    float new_surfaceTensionThreshold,
    float new_surfaceTensionCoeff,
    float new_interfaceTensionThreshold,
    float new_interfaceTensionCoeff,
    float new_diffusionCoeff
    ) {
    // Update parameters
    fluid1_density = new_fluid1_density;
    fluid1_viscosity = new_fluid1_viscosity;
    fluid2_density = new_fluid2_density;
    fluid2_viscosity = new_fluid2_viscosity;
    gravity = new_gravity;
    fluid1_idealGasConstant = new_fluid1_idealGasConstant;
    fluid2_idealGasConstant = new_fluid2_idealGasConstant;
    surfaceTensionThreshold = new_surfaceTensionThreshold;
    surfaceTensionCoeff = new_surfaceTensionCoeff;
    interfaceTensionThreshold = new_interfaceTensionThreshold;
    interfaceTensionCoeff = new_interfaceTensionCoeff;
    diffusionCoeff = new_diffusionCoeff;

    // Only update h and kernel coefficients if h has changed
    if (h != new_h) {
        h = new_h;
        updateKernelCoefficients();
    }

    // Update rest density, viscosity

    for (int i = 0; i < m_fluids[0]->numParticles; i++) {

        m_particles[i]->restDensity = fluid1_density;

    }

    for (int i = 0; i < m_fluids[1]->numParticles; i++) {

        m_particles[m_fluids[0]->numParticles + i]->restDensity = fluid2_density;

    }

    m_fluids[0]->viscosity = fluid1_viscosity;
    m_fluids[1]->viscosity = fluid2_viscosity;

    // Update the density and pressure of all particles
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(i, m_particles[i]->density, m_particles[i]->restDensity);
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

    // Reinitialize simulation
    init();
}
