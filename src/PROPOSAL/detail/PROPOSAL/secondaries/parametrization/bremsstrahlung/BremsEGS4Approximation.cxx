#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsEGS4Approximation.h"
#include "PROPOSAL/Constants.h"
#include <cmath>

using namespace PROPOSAL;

std::pair<Cartesian3D, Cartesian3D>
secondaries::BremsEGS4Approximation::CalculateDirections(
        const Vector3D& init_direction, double energy, double, const Component&,
        std::vector<double>& rnd){
    // Approximation for bremsstrahlung photons used in EGS4
    auto cosphi = std::cos(primary_lepton_mass / energy);
    auto theta = rnd[0] * 2. * PI;
    auto dir_photon = Cartesian3D(init_direction);
    dir_photon.deflect(cosphi, theta);
    return std::pair<Cartesian3D, Cartesian3D>(init_direction, dir_photon);
}