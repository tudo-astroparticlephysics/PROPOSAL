#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"

namespace PROPOSAL {

struct ContRand {
    ContRand() = default;
    virtual ~ContRand() = default;
    static Interpolant1DBuilder::Definition interpol_def;
    virtual double EnergyRandomize(double, double, double) = 0;
};

} // namespace PROPOSAL
