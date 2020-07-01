#include "PROPOSAL/secondaries/SecondariesCalculator.h"

using namespace PROPOSAL;

using std::get;

vector<Loss::secondary_t> SecondariesCalculator::CalculateSecondaries(
    double primary_energy, Loss::secondary_t loss, const Component& comp,
    vector<double> rnd)
{
    auto type = static_cast<InteractionType>(get<Loss::TYPE>(loss));
    auto it = secondary_generator.find(type);
    if (it != secondary_generator.end())
        return it->second->CalculateSecondaries(
                primary_energy, loss, comp, rnd);
    std::ostringstream s;
    s << "No secondary calculator for interaction type ("
      << Type_Interaction_Name_Map.find(type)->second << ") available.";
    throw std::logic_error(s.str());
}

size_t SecondariesCalculator::RequiredRandomNumbers(InteractionType type) const noexcept
{
    return secondary_generator.at(type)->RequiredRandomNumbers();
}
