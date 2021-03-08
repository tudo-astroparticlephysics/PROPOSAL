#include "PROPOSAL/scattering/stochastic_deflection/pairProd/KelnerPairProduction.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
#include <random>
#include "PROPOSAL/math/MathMethods.h" 

using namespace PROPOSAL;

DirectionChangeAngular
stochastic_deflection::KelnerPairProduction::CalculateStochasticDeflection(
        double e_i, double e_f, std::vector<double> const& rnd) const
{
    // All energies should be in units of GeV
    e_i = e_i / 1000.0;
    e_f = e_f / 1000.0;
    auto electron_mass = ME / 1000.0;
    auto muon_mass = mass / 1000.0;

    // Muon values
    auto n = -1.0;
    auto a = 8.9 / 10000.0; 
    auto b = 1.5 / 100000.0;
    auto c = 0.032;
    auto d = 1.0;
    auto e = 0.1;

    auto nu = (e_i - e_f) / (e_i - muon_mass);
    auto min = std::min(a * pow(nu, 0.25) * (1.0 + b * e_i) + c * nu / (nu + d), e);
    auto rms_theta = (2.3 + log(e_i)) * pow(1.0 - nu, n) / e_i * pow(nu - 2.0 * electron_mass / e_i, 2) / (nu * nu) * min;
    
    auto lambda = 1 / (rms_theta * rms_theta);
    // Need sqrt because of sampling in theta^2
    auto theta_muon = std::sqrt(SampleFromExponential(rnd[0], lambda));

    return DirectionChangeAngular { theta_muon, 2 * PI * rnd[1] };
}