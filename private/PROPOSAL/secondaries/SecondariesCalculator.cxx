#include "PROPOSAL/secondaries/SecondariesCalculator.h"

using namespace PROPOSAL;

using std::get;

vector<Loss::secondary_t> SecondariesCalculator::CalculateSecondaries(
    double primary_energy, Loss::secondary_t loss, const Component& comp,
    vector<double> rnd)
{
    auto type = static_cast<InteractionType>(get<Loss::TYPE>(loss));
    return secondary_generator[type]->CalculateSecondaries(
        primary_energy, loss, comp, rnd);
}

size_t SecondariesCalculator::RequiredRandomNumbers(InteractionType type) const noexcept
{
    return secondary_generator.at(type)->RequiredRandomNumbers();
}
