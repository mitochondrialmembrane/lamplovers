# Fluid Simulation: Particle-Based Fluid-Fluid Interaction

This project implements a particle-based fluid-fluid interaction simulation in C++, based on the paper [*"Particle-Based Fluid-Fluid Interaction"*](https://dl.acm.org/doi/pdf/10.1145/1073368.1073402) by Müller et al. The simulation leverages **Smoothed Particle Hydrodynamics (SPH)** to enable mutual interaction between different fluid types (e.g., water and oil). Results are exported as Alembic files for high-quality rendering.

## Features

- Real-time particle-based fluid simulation
- Support for multiple interacting fluid types
- Efficient SPH computation using the **Eigen** library
- Interactive GUI using **Qt**
- Export to `.abc` format with the **Alembic** library for rendering in Blender or other 3D software

## Results:

https://github.com/user-attachments/assets/7a834913-e9db-474e-aba1-b7aa2fa63d50

## Getting Started

### Prerequisites

Make sure you have the following installed:

- C++17 compiler
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Qt 5 or 6](https://www.qt.io/)
- [Alembic](https://github.com/alembic/alembic)

### Build Instructions

```bash
# Clone the repository
git clone https://github.com/yourusername/fluid-simulation.git
cd fluid-simulation

# Create and enter the build directory
mkdir build && cd build

# Configure the project with CMake
cmake ..

# Build the application
make
```
### Usage

Use the GUI to control simulation parameters like:

- **Fluid types and interaction properties**
- **Kernel radius**
- **Time step and iteration count**

Other features:

- Click the **Export** button to generate an Alembic (`.abc`) file.
- Import the `.abc` file into Blender for final rendering.
