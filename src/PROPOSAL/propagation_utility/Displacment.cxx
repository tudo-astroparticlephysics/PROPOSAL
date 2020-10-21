#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

std::unique_ptr<Interpolant1DBuilder::Definition> Displacement::interpol_def = nullptr;

double Displacement::FunctionToIntegral(double energy)
{
    auto result = 0.0;
    for (auto& cr : cross_list)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}
