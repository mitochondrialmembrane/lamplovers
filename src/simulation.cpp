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
    h(settings.value("Parameters/smoothingLength", 0.15f).toFloat()),
    surfaceTensionThreshold(settings.value("Parameters/surfaceTensionThreshold", 0.5f).toFloat()),
    surfaceTensionCoeff(settings.value("Parameters/surfaceTensionCoeff", 20.0f).toFloat()),
    interfaceTensionThreshold(settings.value("Parameters/interfaceTensionThreshold", 0.5f).toFloat()),
    interfaceTensionCoeff(settings.value("Parameters/interfaceTensionCoeff", 20.0f).toFloat()),
    diffusionCoeff(settings.value("Parameters/diffusionCoeff", 0.0001f).toFloat()),
    gravity(Vector3d(settings.value("Parameters/gravity_x", 0.0f).toFloat(),
                     settings.value("Parameters/gravity_y", -10.0f).toFloat(),
                     settings.value("Parameters/gravity_z", 0.0f).toFloat())),
    radius(0.4),
    ceiling(1),
    coneTop(10000)
{
    // Initialize kernel coefficients
    updateKernelCoefficients();
    m_exporter.init("simulation.abc", 100);
    
    // Load fluid configurations from settings
    initFluidsFromSettings();

    // Initialize the fluid indices array with zeros
    m_fluidStartIndices.resize(m_fluids.size() + 1, 0);
}

Simulation::~Simulation() {
    // Clean up particles
    for (auto particle : m_particles) {
        delete particle;
    }
    m_particles.clear();
}

void Simulation::initFluidsFromSettings() {
    // Get the number of fluids from settings
    int numFluids = _settings.value("Parameters/num_fluids", 2).toInt();
    
    // Loop through each fluid and load its parameters
    for (int i = 0; i < numFluids; i++) {
        std::unique_ptr<Fluid> fluid = std::make_unique<Fluid>();
        
        QString prefix = QString("Fluid%1/").arg(i+1);
        
        // Load fluid properties
        fluid->name = _settings.value(prefix + "name", QString("Fluid %1").arg(i+1)).toString().toStdString();
        fluid->density = _settings.value(prefix + "density", 1000.0f).toFloat();
        fluid->restDensity = fluid->density; // Initially same as density
        fluid->viscosity = _settings.value(prefix + "viscosity", 100.0f).toFloat();
        fluid->mass = _settings.value(prefix + "mass", 1.0f).toFloat();
        fluid->gasConstant = _settings.value(prefix + "idealGasConstant", 40.0f).toFloat();
        fluid->colorI = _settings.value(prefix + "colorI", 0.0f).toFloat();
        fluid->colorS = _settings.value(prefix + "colorS", 1.0f).toFloat();
        fluid->numParticles = 0; // Will be set during initialization
        fluid->temperatureDependent = _settings.value(prefix + "temperatureDependent", false).toBool();
        fluid->temperatureConstant = _settings.value(prefix + "temperatureConstant", 10000.0f).toFloat();
        
        // Add the fluid to the simulation
        m_fluids.push_back(std::move(fluid));
        
        // Create a point cloud for this fluid
        m_pointclouds.push_back(std::make_unique<PointCloud>(20.0f)); // Default point size
    }
}

void Simulation::updateFluidIndices() {
    // Ensure the indices array is the correct size (fluids + 1 for the end index)
    m_fluidStartIndices.resize(m_fluids.size() + 1, 0);
    
    // Calculate starting indices for each fluid for faster access
    m_fluidStartIndices[0] = 0; // First fluid starts at index 0
    for (size_t i = 0; i < m_fluids.size(); i++) {
        m_fluidStartIndices[i+1] = m_fluidStartIndices[i] + m_fluids[i]->numParticles;
    }
    
    // Verify that the last index matches the total particle count
    if (m_fluidStartIndices.back() != static_cast<int>(m_particles.size())) {
        std::cerr << "Warning: Fluid indices don't match particle count. Expected " 
                  << m_particles.size() << " but got " << m_fluidStartIndices.back() << std::endl;
    }
}

void Simulation::init() {
    // Initialize fluids and their particles
    if (m_fluids.empty()) {
        std::cerr << "Error: No fluids defined in settings!" << std::endl;
        return;
    }
    
    // Create two zones for the fluids (can be extended for more complex scenarios)
    int d = 12; // Density of particles
    
    float totalHeight = ceiling * d;
    float heightPerFluid = totalHeight / m_fluids.size();
    
    // Generate particles for each fluid
    for (size_t fluidIndex = 0; fluidIndex < m_fluids.size(); fluidIndex++) {
        float startHeight = fluidIndex * heightPerFluid / d;
        float endHeight = (fluidIndex + 1) * heightPerFluid / d;
        
        std::vector<Point> fluidPoints;
        
        #pragma omp parallel
        {
            std::vector<Vector3d> localPositions;
            
            #pragma omp for collapse(3)
            for (int i = -radius * d; i < radius * d; i++) {
                for (int j = -radius * d; j < radius * d; j++) {
                    for (int k = startHeight * d; k < endHeight * d; k++) {
                        Vector3d pos((float) i / d, (float) k / d, (float) j / d);
                        if (checkCollision(pos) == Vector3d(0,0,0)) {
                            localPositions.push_back(pos);
                        }
                    }
                }
            }
            
            // Add particles to the main list
            #pragma omp critical
            {
                for (const auto& pos : localPositions) {
                    Particle* newParticle = new Particle{
                        pos, 
                        Vector3d(0, 0, 0), // Initial velocity
                        0, // Initial pressure
                        m_fluids[fluidIndex]->density, 
                        m_fluids[fluidIndex]->restDensity, 
                        10, // Initial temperature
                        m_fluids[fluidIndex].get()
                    };
                    
                    m_particles.push_back(newParticle);
                    fluidPoints.push_back(Point{pos, 10}); // Initial temperature
                }
            }
        }
        
        // Set the number of particles for this fluid
        m_fluids[fluidIndex]->numParticles = fluidPoints.size();
        
        // Initialize the point cloud for this fluid
        m_pointclouds[fluidIndex]->init(fluidPoints, fluidIndex % 2); // Alternate colors for now
    }
    
    // Calculate and store the starting indices for each fluid
    updateFluidIndices();
    
    initGround();
    initBox();
}

void Simulation::update(double seconds) {
    seconds = std::min(seconds, 0.01);  // Cap maximum time step
    
    // Safety check - don't proceed if there are no particles
    if (m_particles.empty()) return;

    std::vector<Vector3d> accelerations;
    accelerations.resize(m_particles.size());

    // Calculate initial accelerations
    #pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        if (!m_particles[i]) continue; // Skip null particles
        accelerations[i] = (fPressure(i) + fViscosity(i)) / m_particles[i]->density + gravity;
    }

    std::vector<Vector3d> newAccelerations = accelerations;
    
    // Update positions using Leapfrog integration
    #pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->position += m_particles[i]->velocity * seconds + 0.5 * accelerations[i] * seconds * seconds;
    }

    // Update densities and pressures
    #pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(i, m_particles[i]->density, m_particles[i]->restDensity);
    }

    // Calculate new accelerations with all forces
    #pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        newAccelerations[i] = (fPressure(i) + fViscosity(i) + fSurfaceTension(i) + fInterfaceTension(i)) / m_particles[i]->density + gravity;
    }

    // Update velocities and handle collisions
    #pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->velocity += 0.5 * (accelerations[i] + newAccelerations[i]) * seconds;
        evaluateCollisions(i);
    }

    // Final density and pressure update
    #pragma omp parallel for
    for (int i = 0; i < m_particles.size(); i++) {
        m_particles[i]->density = density_S(i);
        m_particles[i]->pressure = calculatePressure(i, m_particles[i]->density, m_particles[i]->restDensity);
    }

    // Temperature diffusion and effects
    
    // Heat sources and sinks (bottom and top)
    for (int i = 0; i < m_particles.size(); i++) {
        if (m_particles[i]->position[1] < 0.3) {
            m_particles[i]->temperature = fmin(m_particles[i]->temperature + 0.1, 50);
        }

        if (m_particles[i]->position[1] > 0.7) {
            m_particles[i]->temperature = fmax(m_particles[i]->temperature - 0.1, 1);
        }
    }

    // Heat diffusion
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
        m_particles[i]->temperature = fmax(1.0, fmin(m_particles[i]->temperature, 50.0));
        
        // Update rest density for temperature-dependent fluids
        if (m_particles[i]->fluid->temperatureDependent) {
            m_particles[i]->restDensity = m_particles[i]->fluid->temperatureConstant / m_particles[i]->temperature;
        }
    }

    // Update point clouds for visualization using fluid indices for more efficient access
    std::vector<std::vector<Point>> fluidPoints(m_fluids.size());
    std::vector<std::vector<Vector3d>> fluidPositions(m_fluids.size());
    
    // Pre-allocate space for points to avoid constant reallocation
    for (size_t fluidIdx = 0; fluidIdx < m_fluids.size(); fluidIdx++) {
        fluidPoints[fluidIdx].reserve(m_fluids[fluidIdx]->numParticles);
        fluidPositions[fluidIdx].reserve(m_fluids[fluidIdx]->numParticles);
    }
    
    // Use fluid indices to directly access each fluid's particles
    for (size_t fluidIdx = 0; fluidIdx < m_fluids.size(); fluidIdx++) {
        for (int i = m_fluidStartIndices[fluidIdx]; i < m_fluidStartIndices[fluidIdx+1]; i++) {
            fluidPoints[fluidIdx].push_back(Point{
                m_particles[i]->position, 
                m_particles[i]->temperature
            });
            fluidPositions[fluidIdx].push_back(m_particles[i]->position);
        }
        
        // Update point cloud for rendering
        m_pointclouds[fluidIdx]->setPoints(fluidPoints[fluidIdx]);
    }
    
    // Update exporter
    if (m_fluids.size() >= 2) {  // Ensure backward compatibility with exporter
        m_exporter.addFrame(fluidPositions[0], fluidPositions[1]);
    }
}

void Simulation::draw(Shader *shader) {
    m_shape.draw(shader);
    m_ground.draw(shader);
    
    // Draw all fluid point clouds
    for (auto& pointcloud : m_pointclouds) {
        pointcloud->draw(shader);
    }
    
    m_box.draw(shader);
}

void Simulation::initGround() {
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
    // Using reduction to safely accumulate the sum in parallel
    #pragma omp parallel for reduction(+:out)
    for (int j = 0; j < m_particles.size(); j++) {
        double r_squared = (m_particles[i]->position - m_particles[j]->position).squaredNorm() + 0.00001;
        if (r_squared < h * h) {
            out += m_particles[j]->fluid->mass * Wpoly6Coeff * pow((h * h - r_squared), 3);
        }
    }

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
            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * 
                                 (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));
            out += diffusionCoeff * (particleJ->fluid->mass / particleJ->density * W_Laplacian) * 
                                    (particleJ->temperature - particleI->temperature);
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
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-WspikyGradCoeff * pow((h - r_norm), 2) / r_norm) * r;
            out += particleJ->fluid->mass * (particleJ->pressure + pressureI) / 
                  (2 * particleJ->density) * W_Grad;
        }
    }
    return -out;
}

Vector3d Simulation::fViscosity(int i) {
    Particle *particleI = m_particles[i];
    Vector3d out(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        Particle *particleJ = m_particles[j];
        Vector3d r = particleI->position - particleJ->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            double W_Laplacian = WviscosityLaplacianCoeff * (h - r_norm);
            double avgViscosity = (particleJ->fluid->viscosity + particleI->fluid->viscosity) / 2.0;
            out += avgViscosity * (particleJ->fluid->mass / particleJ->density * W_Laplacian) * 
                  (particleJ->velocity - particleI->velocity);
        }
    }

    return out;
}

Vector3d Simulation::fSurfaceTension(int i) {
    Vector3d normal(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-Wpoly6GradCoeff * pow((h - r_norm), 2)) * r;
            normal += (particleJ->fluid->mass / particleJ->density) * W_Grad;
        }
    }
    
    if (normal.norm() < surfaceTensionThreshold) return Vector3d(0,0,0);

    double C_S_Laplacian = 0;
    for (int j = 0; j < m_particles.size(); j++) {
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_squaredNorm = r.squaredNorm() + 0.00001;
        double h_squared = h * h;
        if (r_squaredNorm < h_squared) {
            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * 
                                (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));
            C_S_Laplacian += (particleJ->fluid->mass / particleJ->density) * W_Laplacian;
        }
    }
    return (-surfaceTensionCoeff * C_S_Laplacian) * normal.normalized();
}

Vector3d Simulation::fInterfaceTension(int i) {
    Vector3d normal(0,0,0);
    for (int j = 0; j < m_particles.size(); j++) {
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_norm = r.norm() + 0.00001;
        if (r_norm < h) {
            Vector3d W_Grad = (-Wpoly6GradCoeff * pow((h - r_norm), 2)) * r;
            normal += particleJ->fluid->colorI * (particleJ->fluid->mass / particleJ->density) * W_Grad;
        }
    }
    
    if (normal.norm() < interfaceTensionThreshold) return Vector3d(0,0,0);

    double C_I_Laplacian = 0;
    for (int j = 0; j < m_particles.size(); j++) {
        Particle *particleJ = m_particles[j];
        Vector3d r = m_particles[i]->position - m_particles[j]->position;
        double r_squaredNorm = r.squaredNorm() + 0.00001;
        double h_squared = h * h;
        if (r_squaredNorm < h_squared) {
            double W_Laplacian = Wpoly6LaplacianCoeff * (h_squared - r_squaredNorm) * 
                                (r_squaredNorm - 0.75 * (h_squared - r_squaredNorm));
            C_I_Laplacian += particleJ->fluid->colorI * (particleJ->fluid->mass / particleJ->density) * W_Laplacian;
        }
    }
    return (-interfaceTensionCoeff * C_I_Laplacian) * normal.normalized();
}

void Simulation::evaluateCollisions(int i) {
    Vector3d displacement = checkCollision(m_particles[i]->position);

    if (displacement.norm() > 0.0001) {
        m_particles[i]->position += displacement;

        // Calculate damped reflection with friction
        float dampingCoeff = 0.7; // Lower than 1.0 to dissipate energy
        Vector3d normal = displacement.normalized();
        float vn = m_particles[i]->velocity.dot(normal);

        if (vn < 0) { // Only dampen if moving into the boundary
            Vector3d vTangent = m_particles[i]->velocity - vn * normal;
            float friction = 0.1; // Tangential friction coefficient

            m_particles[i]->velocity = vTangent * (1.0f - friction) - dampingCoeff * vn * normal;
        }
    }
}

Vector3d Simulation::checkCollision(Vector3d pos) {
    float r = sqrt(pos[0] * pos[0] + pos[2] * pos[2]);
    float radiusAtY = (coneTop - pos[1]) / coneTop * radius;
    Vector3d returnVector(0,0,0);
    if (pos[1] < 0) returnVector += Vector3d(0,1,0) * -pos[1];
    if (pos[1] > ceiling) returnVector += Vector3d(0,-1,0) * (pos[1] - ceiling);
    if (r > radiusAtY) returnVector += (Vector3d(-pos[0], 0, -pos[2]).normalized() * coneTop - Vector3d(0, -radius, 0)).normalized() * (r - radiusAtY);
    return returnVector;
}

void Simulation::updateParameters(const std::map<QString, float>& paramValues) {
    // Update global parameters
    for (const auto& [param, value] : paramValues) {
        if (param == "smoothingLength") {
            if (h != value) {
                h = value;
                updateKernelCoefficients();
            }
        } else if (param == "gravity_x") {
            gravity[0] = value;
        } else if (param == "gravity_y") {
            gravity[1] = value;
        } else if (param == "gravity_z") {
            gravity[2] = value;
        } else if (param == "surfaceTensionThreshold") {
            surfaceTensionThreshold = value;
        } else if (param == "surfaceTensionCoeff") {
            surfaceTensionCoeff = value;
        } else if (param == "interfaceTensionThreshold") {
            interfaceTensionThreshold = value;
        } else if (param == "interfaceTensionCoeff") {
            interfaceTensionCoeff = value;
        } else if (param == "diffusionCoeff") {
            diffusionCoeff = value;
        }
    }
    
    // Update fluid-specific parameters
    for (size_t i = 0; i < m_fluids.size(); i++) {
        QString prefix = QString("fluid%1_").arg(i+1);
        
        for (const auto& [param, value] : paramValues) {
            if (param.startsWith(prefix)) {
                QString fluidParam = param.mid(prefix.length());
                
                if (fluidParam == "density") {
                    m_fluids[i]->density = value;
                    
                    // Update rest density for non-temperature dependent particles
                    if (!m_fluids[i]->temperatureDependent) {
                        m_fluids[i]->restDensity = value;
                        
                        // Update particle rest densities using fluid indices
                        for (int p = m_fluidStartIndices[i]; p < m_fluidStartIndices[i+1]; p++) {
                            m_particles[p]->restDensity = value;
                        }
                    }
                } else if (fluidParam == "viscosity") {
                    m_fluids[i]->viscosity = value;
                } else if (fluidParam == "idealGasConstant") {
                    m_fluids[i]->gasConstant = value;
                } else if (fluidParam == "colorI") {
                    m_fluids[i]->colorI = value;
                }
            }
        }
    }
    
    // Recalculate all densities and pressures based on new parameters
    #pragma omp parallel for
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

void Simulation::reinitialize() {
    // Clean up existing particles
    for (auto particle : m_particles) {
        delete particle;
    }
    m_particles.clear();
    
    // Reset particle counters in fluids
    for (auto& fluid : m_fluids) {
        fluid->numParticles = 0;
    }
    
    // Reset fluid indices
    m_fluidStartIndices.clear();
    m_fluidStartIndices.resize(m_fluids.size() + 1, 0);
    
    // Reinitialize simulation
    init();
}

void Simulation::createFluidParticles(Fluid* fluid, const Eigen::Vector3d& startPos, 
                                     const Eigen::Vector3d& dimensions, double spacing) {
    std::vector<Vector3d> positions;
    std::vector<Point> points;
    
    // Calculate number of particles in each dimension
    int nx = static_cast<int>(dimensions.x() / spacing);
    int ny = static_cast<int>(dimensions.y() / spacing);
    int nz = static_cast<int>(dimensions.z() / spacing);
    
    // Create a grid of particles
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Vector3d pos = startPos + Vector3d(
                    i * spacing, 
                    j * spacing, 
                    k * spacing
                );
                
                // Check if position is valid (inside boundaries)
                if (checkCollision(pos) == Vector3d(0,0,0)) {
                    positions.push_back(pos);
                    points.push_back(Point{pos, 10}); // Initial temperature
                }
            }
        }
    }
    
    // Find the fluid index
    int fluidIndex = -1;
    for (size_t i = 0; i < m_fluids.size(); i++) {
        if (m_fluids[i].get() == fluid) {
            fluidIndex = i;
            break;
        }
    }
    
    if (fluidIndex == -1) {
        std::cerr << "Error: Fluid not found in simulation" << std::endl;
        return;
    }
    
    // Create particles and add them to the simulation at the correct position
    // (after the last particle of the previous fluid)
    int insertPosition = m_fluidStartIndices[fluidIndex+1];
    
    // Create temporary storage for new particles
    std::vector<Particle*> newParticles;
    newParticles.reserve(positions.size());
    
    for (const auto& pos : positions) {
        Particle* newParticle = new Particle{
            pos, 
            Vector3d(0, 0, 0), // Initial velocity
            0, // Initial pressure
            fluid->density, 
            fluid->restDensity, 
            10, // Initial temperature
            fluid
        };
        
        newParticles.push_back(newParticle);
    }
    
    // Insert the new particles at the correct position
    m_particles.insert(m_particles.begin() + insertPosition, newParticles.begin(), newParticles.end());
    
    // Update fluid's particle count
    fluid->numParticles += positions.size();
    
    // Update start indices for all subsequent fluids
    for (size_t i = fluidIndex + 1; i <= m_fluids.size(); i++) {
        m_fluidStartIndices[i] += positions.size();
    }
    
    // Find the pointcloud for this fluid and update it
    if (fluidIndex >= 0 && fluidIndex < m_pointclouds.size()) {
        // Get all particles for this fluid
        std::vector<Point> allPoints;
        allPoints.reserve(fluid->numParticles);
        
        for (int i = m_fluidStartIndices[fluidIndex]; i < m_fluidStartIndices[fluidIndex+1]; i++) {
            allPoints.push_back(Point{m_particles[i]->position, m_particles[i]->temperature});
        }
        
        m_pointclouds[fluidIndex]->setPoints(allPoints);
    }
}

double Simulation::calculateTimeStep() {
    double maxVelocity = 0.0;
    double minDistance = h; // Start with smoothing length
    
    // Find maximum velocity and minimum particle distance
    for (int i = 0; i < m_particles.size(); i++) {
        double velMag = m_particles[i]->velocity.norm();
        maxVelocity = std::max(maxVelocity, velMag);
        
        // Find nearest neighbor distance (can be optimized with spatial hash)
        for (int j = i+1; j < m_particles.size(); j++) {
            double dist = (m_particles[i]->position - m_particles[j]->position).norm();
            if (dist > 0.0001) minDistance = std::min(minDistance, dist);
        }
    }
    
    // CFL condition with safety factor
    double safetyFactor = 0.4;
    double timeStep = safetyFactor * minDistance / (maxVelocity + 1e-6);
    
    // Clamp to reasonable range
    return std::min(std::max(timeStep, 0.001), 0.01);
}

void Simulation::addFluid(const Fluid& fluidProperties) {
    // Create a new fluid
    std::unique_ptr<Fluid> newFluid = std::make_unique<Fluid>(fluidProperties);
    
    // Add it to the simulation
    m_fluids.push_back(std::move(newFluid));
    
    // Create a point cloud for this fluid
    m_pointclouds.push_back(std::make_unique<PointCloud>(20.0f));
    
    // Update fluid indices array
    m_fluidStartIndices.push_back(m_fluidStartIndices.back());
}
