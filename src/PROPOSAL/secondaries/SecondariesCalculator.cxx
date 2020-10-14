#include "PROPOSAL/secondaries/SecondariesCalculator.h"

using namespace PROPOSAL;

std::vector<Loss::secondary_t> SecondariesCalculator::CalculateSecondaries(
        double initial_energy, Loss::secondary_t loss, const Component& comp,
        vector<double> rnd)
{
    auto type = static_cast<InteractionType>(std::get<Loss::TYPE>(loss));
    auto it = secondary_generator.find(type);
    if (it != secondary_generator.end())
        return it->second->CalculateSecondaries(
                initial_energy, loss, comp, rnd);
    std::ostringstream s;
    s << "No secondary calculator for interaction type ("
        << Type_Interaction_Name_Map.find(type)->second << ") available.";
    throw std::logic_error(s.str());
}
