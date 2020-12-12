#include "PROPOSAL/secondaries/SecondariesCalculator.h"

using namespace PROPOSAL;

std::vector<ParticleState> SecondariesCalculator::CalculateSecondaries(
        StochasticLoss loss, const Component& comp, std::vector<double> &rnd)
{
    auto type = static_cast<InteractionType>(loss.type);
    auto it = secondary_generator.find(type);
    if (it != secondary_generator.end())
        return it->second->CalculateSecondaries(loss, comp, rnd);
    std::ostringstream s;
    s << "No secondary calculator for interaction type ("
        << Type_Interaction_Name_Map.find(type)->second << ") available.";
    throw std::logic_error(s.str());
}
