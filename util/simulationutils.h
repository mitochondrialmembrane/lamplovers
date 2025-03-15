
#include <Eigen/Dense>

using namespace Eigen;

namespace SimulationUtils
{
double tetrahedronVolume(const Vector3d& v1, const Vector3d& v2,
                         const Vector3d& v3, const Vector3d& v4);
}
