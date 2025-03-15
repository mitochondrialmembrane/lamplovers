#include "simulationutils.h"

#include <iostream>
#include <vector>

using namespace Eigen;

namespace SimulationUtils
{
double tetrahedronVolume(const Vector3d& v1, const Vector3d& v2,
                         const Vector3d& v3, const Vector3d& v4) {
    Eigen::Vector3d v2_v1 = v2 - v1;
    Eigen::Vector3d v3_v1 = v3 - v1;
    Eigen::Vector3d v4_v1 = v4 - v1;

    // Construct the 3x3 matrix with difference vectors as columns
    Eigen::Matrix3d M;
    M << v2_v1, v3_v1, v4_v1;

    // Compute volume using determinant formula
    return std::abs(M.determinant()) / 6.0;
}

}
