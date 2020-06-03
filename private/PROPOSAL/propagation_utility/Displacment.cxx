#include <algorithm>

#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;
using std::max;

template <typename Cross>
double Displacement::FunctionToIntegral(Cross&& cross, double energy)
{
    auto result = 0.0;
    for (auto& cr : cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}

Interpolant1DBuilder::Definition Displacement::interpol_def{};

namespace PROPOSAL {
Interpolant1DBuilder::Definition displacement_interpol_def;
} // namespace PROPOSAL
