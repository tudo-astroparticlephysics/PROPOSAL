#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/TsaiApproximationBremsstrahlung.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

DirectionChangeAngular 
stochastic_deflection::TsaiApproximationBremsstrahlung::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd) const
{
    auto epsilon = e_i - e_f;
    auto theta_star = 1.0;
    auto r_max = std::min(1.0, e_f/epsilon) * e_i * theta_star / mass;
    auto a = rnd[0] * r_max*r_max / (1.0 + r_max*r_max);
    auto r = std::sqrt(a / (1.0-a));
    auto theta_photon = mass / e_i * r;
    auto theta_muon = epsilon / e_f * theta_photon;

    return  DirectionChangeAngular {theta_muon, 2 * PI * rnd[1]};
}