
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

Displacement::Displacement(const CrossSectionList& cross)
    : cross(cross)
{
    if (cross.size() < 1)
        throw std::invalid_argument("At least one crosssection is required.");
}

double Displacement::FunctionToIntegral(double energy)
{
    auto result = 0.0;
    for (const auto& cr : cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}
