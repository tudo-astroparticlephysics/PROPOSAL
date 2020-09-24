#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

Interpolant1DBuilder::Definition Displacement::interpol_def(1000);

double Displacement::FunctionToIntegral(double energy)
{
    auto result = 0.0;
    for (auto& cr : cross_list)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}
