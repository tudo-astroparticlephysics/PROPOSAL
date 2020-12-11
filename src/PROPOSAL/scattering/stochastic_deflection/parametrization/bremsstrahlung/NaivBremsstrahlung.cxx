#include "PROPOSAL/scattering/stochastic_deflection/parametrization/bremsstrahlung/NaivBremsstrahlung.h"

using namespace PROPOSAL;

std::array<double, 2>
stochastic_deflection::NaivBremsstrahlung::CalculateDeflection(
    StochasticLoss loss, const Component&, std::vector<double>&)
{
    return std::array<double, 2> { 1.23, 3.45 };
}
