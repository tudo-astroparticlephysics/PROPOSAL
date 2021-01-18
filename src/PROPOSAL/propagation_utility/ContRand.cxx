#include "PROPOSAL/propagation_utility/ContRand.h"

using namespace PROPOSAL;

Interpolant1DBuilder::Definition ContRand::interpol_def = { 200 };

double ContRand::FunctionToIntegral(double energy)
{
    assert(energy >= 0);
    double sum = 0.0;
    for (auto& crosssections : cross_list)
        sum += crosssections->CalculatedE2dx(energy);
    return disp->FunctionToIntegral(energy) * sum;
}
