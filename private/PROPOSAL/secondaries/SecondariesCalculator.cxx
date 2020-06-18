#include "PROPOSAL/secondaries/SecondariesCalculator.h"

using namespace PROPOSAL;

using std::get;

vector<Loss::secondary_t> SecondariesCalculator::CalculateSecondaries(
    double primary_energy, Loss::secondary_t loss, const Component& comp,
    vector<double> rnd)
{
    auto& calculator = secondary_generator[static_cast<int>(get<Loss::TYPE>(loss))];
    return calculator->CalculateSecondaries(primary_energy, loss, comp, rnd);
}

size_t SecondariesCalculator::requiredRandomNumbers(int type)
{
    return secondary_generator[type]->RequiredRandomNumbers();
}
