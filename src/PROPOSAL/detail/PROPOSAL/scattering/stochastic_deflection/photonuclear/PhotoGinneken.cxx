#include "PROPOSAL/scattering/stochastic_deflection/photonuclear/PhotoGinneken.h"
#include "PROPOSAL/Constants.h"
#include <cmath> 
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

UnitSphericalVector 
stochastic_deflection::PhotoGinneken::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd, size_t) const 
{
    // All energies should be in units of GeV
    e_i = e_i / 1000.0;
    e_f = e_f / 1000.0;
    auto muon_mass = mass / 1000.0;

    auto nu = (e_i - e_f) / (e_i - muon_mass);
    auto rms_theta = (0.39 / (e_i * (1 - nu))) * pow(std::sqrt(e_i) * nu *(1 - nu), 0.17) * (1 - 0.135 / (e_i * nu));

    auto lambda = 1 / (rms_theta * rms_theta);
    // Need sqrt because of sampling in theta^2
    auto theta_muon = std::sqrt(SampleFromExponential(rnd[0], lambda));

    return UnitSphericalVector { theta_muon, 2 * PI * rnd[1] };
}