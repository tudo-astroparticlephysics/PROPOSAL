#include "PROPOSAL/scattering/stochastic_deflection/StochasticDeflectionCalculator.h"

using namespace PROPOSAL;

std::array<double, 2> StochasticDeflectionCalculator::CalculateDeflection(
        StochasticLoss loss, const Component& comp, std::vector<double> &rnd)
{
    auto type = static_cast<InteractionType>(loss.type);
    auto it = m.find(type);
    if (it != m.end())
        return it->second->CalculateStochasticDeflection(loss, comp, rnd);
    std::ostringstream s;
    s << "No stochastic deflection calculator for interaction type ("
        << Type_Interaction_Name_Map.find(type)->second << ") available.";
    throw std::logic_error(s.str());
}
