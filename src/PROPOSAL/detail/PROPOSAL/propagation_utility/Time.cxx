#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

Interpolant1DBuilder::Definition Time::interpol_def = { 1000 };

Time::Time(std::shared_ptr<Displacement> _disp, double _mass)
        : disp(_disp)
        , mass(_mass)
        , hash(0)
{
    hash_combine(hash, mass);
}

double Time::FunctionToIntegral(double energy)
{
    assert(energy > mass);
    auto square_momentum = (energy - mass) * (energy + mass);
    auto particle_momentum = std::sqrt(square_momentum);
    return disp->FunctionToIntegral(energy) * energy
        / (particle_momentum * SPEED);
}
