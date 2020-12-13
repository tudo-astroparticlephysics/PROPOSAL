#include "PROPOSAL/scattering/stochastic_deflection/StochasticDeflectionCalculator.h"

using namespace PROPOSAL;

std::array<double, 2> StochasticDeflectionCalculator::CalculateDeflection(
        StochasticLoss loss, const Component& comp, std::vector<double> &rnd) const
{
    auto it = m.find(static_cast<InteractionType>(loss.type));
    if (it != m.end())
        return it->second->CalculateStochasticDeflection(loss, comp, rnd);
    return std::array<double, 2>{0., 0.};
}
