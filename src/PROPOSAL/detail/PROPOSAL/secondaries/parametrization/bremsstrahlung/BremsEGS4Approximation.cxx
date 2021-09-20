#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsEGS4Approximation.h"
#include "PROPOSAL/Constants.h"
#include <cmath>

using namespace PROPOSAL;

std::pair<Cartesian3D, Cartesian3D>
secondaries::BremsEGS4Approximation::CalculateDirections(
        const Vector3D& init_direction, double energy, double photon_energy,
        const Component&, std::vector<double>& rnd){
    // Approximation for bremsstrahlung photons used in EGS4
    auto cosphi = std::cos(primary_lepton_mass / energy);
    auto theta = rnd[0] * 2. * PI;
    auto dir_photon = Cartesian3D(init_direction);
    dir_photon.deflect(cosphi, theta);

    // calculate new electron direction using momentum conservation, neglecting
    // the momentum transfer to the nucleus (i.e. p_f = p_i - p_photon)
    auto dir_electron = Cartesian3D(init_direction) *
            std::sqrt((energy + primary_lepton_mass) *
            (energy - primary_lepton_mass)) - dir_photon * photon_energy;
    dir_electron.normalize();

    return std::pair<Cartesian3D, Cartesian3D>(dir_electron, dir_photon);
}