
#include "PROPOSAL/propagation_utility/Decay.h"

using namespace PROPOSAL;

Decay::Decay(double lifetime, double mass, double lower_lim)
    : lifetime(lifetime)
    , mass(mass)
    , lower_lim(lower_lim)
{
}

std::unique_ptr<Interpolant1DBuilder::Definition> Decay::interpol_def = nullptr;
