#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"

#include <memory>

namespace PROPOSAL {

struct CrossSectionInterpolantBase {
    static Interpolant2DBuilder::Definition dNdx_def;
    static Interpolant1DBuilder::Definition dEdx_def;
    static Interpolant1DBuilder::Definition dE2dx_def;
};


} // namespace PROPOSAL
