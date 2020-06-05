
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
#include <cassert>

using namespace PROPOSAL;

Decay::Decay(double lifetime, double mass, double lower_lim)
    : lifetime(lifetime)
    , mass(mass)
    , lower_lim(lower_lim)
{
}

Interpolant1DBuilder::Definition Decay::interpol_def{};

template <typename Disp>
double Decay::FunctionToIntegral(Disp&& disp, double energy)
{
    assert(!std::isinf(lifetime));
    assert(energy >= mass);
    double square_momentum = (energy - mass) * (energy + mass);
    double aux = SPEED * std::sqrt(square_momentum) / mass;
    return disp.FunctionToIntegral(energy) / aux;
}
