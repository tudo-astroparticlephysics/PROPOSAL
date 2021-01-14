#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

std::array<double, 2>
stochastic_deflection::NaivBremsstrahlung::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd) const
{

    // -------- Bremsstrahlung
    auto epsilon = e_i - e_f;
    auto muon_mass = 105.658;
    auto theta_star = 1.0;
    auto r_max = std::min(1.0, e_f/epsilon) * e_i * theta_star / muon_mass;
    auto a = rnd[0] * r_max*r_max / (1.0 + r_max*r_max);
    auto r = std::sqrt(a / (1.0-a));
    auto theta_photon = muon_mass / e_i * r;
    auto theta_muon = epsilon / e_f * theta_photon;

    return std::array<double, 2> { theta_photon, 2 * PI * rnd[1]};
}
