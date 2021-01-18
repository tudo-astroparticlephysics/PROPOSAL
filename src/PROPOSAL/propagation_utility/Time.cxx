#include "PROPOSAL/propagation_utility/Time.h"

using namespace PROPOSAL;

Interpolant1DBuilder::Definition Time::interpol_def = { 1000 };

double Time::FunctionToIntegral(double energy)
{
    assert(energy > mass);
    auto square_momentum = (energy - mass) * (energy + mass);
    auto particle_momentum = std::sqrt(square_momentum);
    return disp->FunctionToIntegral(energy) * energy
        / (particle_momentum * SPEED);
}
