#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsKochMotz.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include <cmath>
#include <algorithm>
#include <limits>

using namespace PROPOSAL;

double secondaries::BremsKochMotz::CalculatePhiCandidate(double E_0, double rnd) {
    // equation (5)
    return 1. / E_0 * std::sqrt(rnd / (1. - rnd + 1. / std::pow(PI * E_0, 2)));
}

double secondaries::BremsKochMotz::g(double x, double E, double E_0, double Z) {
    // equation (4), but use g(x) directly from equation (2.152) in the EGS5
    // manual since there seems to be a typo in the original paper
    auto r = E / E_0;
    auto m = std::pow((1. - r) / (2. * E_0 * r), 2) +
            std::pow( std::pow(Z, 1./3.) / (111. * (x + 1.)) , 2);
    return 3. * (1. + r * r) - 2. * r - (4. + std::log(m)) *
            ((1. + r * r) - 4. * x * r / std::pow(x + 1., 2));
}

double secondaries::BremsKochMotz::N(double E, double E_0, double Z) {
    // equation (6)
    return 1. / std::max({g(0., E, E_0, Z),
                          g(1., E, E_0, Z),
                          g(std::pow(PI * E_0, 2.), E, E_0, Z)});
}

std::pair<Cartesian3D, Cartesian3D>
secondaries::BremsKochMotz::CalculateDirections(
        const Vector3D& init_direction, double energy, double photon_energy,
        const Component& c, std::vector<double>& rnd){
    // We do not know in advance how many random numbers we need, so use one
    // passed random number to set the seed for a RNG to keep reproducibility
    RandomGenerator internal_rnd_generator;
    internal_rnd_generator.SetSeed(std::numeric_limits<int>::max() * rnd[0]);

    // Based on "Improved bremsstrahlung photon angular sampling in the EGS4
    // code system" (PIRS-0203), algorithm given in chapter (2.1)
    double E = (energy - photon_energy) / ME; // energies in units of ME
    double E_0 = energy / ME;
    auto Z = c.GetNucCharge();
    double phi, x_hat, g_max;
    do {
        phi = CalculatePhiCandidate(E_0, internal_rnd_generator.RandomDouble());
        x_hat = std::pow(E_0 * phi, 2);
        g_max = N(E, E_0, Z) * g(x_hat, E, E_0, Z);
    } while (internal_rnd_generator.RandomDouble() > g_max);
    auto cosphi = std::cos(phi);
    auto theta = rnd[1] * 2. * PI;
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