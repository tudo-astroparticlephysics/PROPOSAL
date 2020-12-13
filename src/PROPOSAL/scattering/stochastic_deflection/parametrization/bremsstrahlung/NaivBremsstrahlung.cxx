#include "PROPOSAL/scattering/stochastic_deflection/parametrization/bremsstrahlung/NaivBremsstrahlung.h"

using namespace PROPOSAL;

std::array<double, 2>
stochastic_deflection::NaivBremsstrahlung::CalculateStochasticDeflection(
    StochasticLoss const &loss, Component const&, std::vector<double> const&) const
{
    return std::array<double, 2> { 1.23, 3.45 };
}
