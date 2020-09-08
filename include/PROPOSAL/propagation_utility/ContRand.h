#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"

namespace PROPOSAL {

struct ContRand {
    ContRand() = default;
    static Interpolant1DBuilder::Definition contrand_interpol_def;
    virtual double EnergyRandomize(double, double, double) = 0;
};

} // namespace PROPOSAL
