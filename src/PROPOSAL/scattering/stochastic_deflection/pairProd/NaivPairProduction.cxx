#include "PROPOSAL/scattering/stochastic_deflection/pairProd/NaivPairProduction.h"
#include "PROPOSAL/Constants.h"
#include "math.h"
#include "random"
#include "PROPOSAL/math/MathMethods.h" 

using namespace PROPOSAL;

std::array<double, 2>
stochastic_deflection::NaivPairProduction::CalculateStochasticDeflection(
        double e_i, double e_f, std::vector<double> const& rnd) const
{
    // ------ PairProd -----------
    auto electron_mass = 0.511; // Proton mass in MeV
    auto muon_mass = 105.558;
    // All energies should be in units of GeV
    e_i = e_i / 1000.0;
    e_f = e_f / 1000.0;
    electron_mass = electron_mass / 1000.0;
    muon_mass = muon_mass / 1000.0;

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

    //TODO: PropagationUtility will call cos() on the first return value
    return std::array<double, 2> { theta_muon, 2 * PI * rnd[1] };
}