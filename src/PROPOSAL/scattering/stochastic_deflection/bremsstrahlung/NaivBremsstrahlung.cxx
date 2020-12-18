#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/NaivBremsstrahlung.h"
#include <iostream>
using namespace PROPOSAL;

std::array<double, 2>
stochastic_deflection::NaivBremsstrahlung::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd) const
{
    auto sign1 = rnd[0] > 0.5 ? 1.: -1.;
    auto sign2 = rnd[1] > 0.5 ? 1.: -1.;
    std::cout << sign1 << std::endl;
    std::cout << sign2 << std::endl;
    return std::array<double, 2> { sign1 * (e_f / e_i) * 0.1, sign2 * (e_f / e_i) * 0.1 };
}
