#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsNoDeflection.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

std::pair<Cartesian3D, Cartesian3D>
        secondaries::BremsNoDeflection::CalculateDirections(
                const Vector3D& init_direction, double, double,
                const Component&, std::vector<double>&) {
    return std::pair<Cartesian3D, Cartesian3D>(init_direction, init_direction);
}