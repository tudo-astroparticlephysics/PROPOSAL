#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/Constants.h"
#include <iostream>
using namespace PROPOSAL;

std::array<double, 2>
stochastic_deflection::NaivBremsstrahlung::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd) const
{
    return std::array<double, 2> { e_f / e_i * 0.1, 2 * PI * rnd[0] };
}
