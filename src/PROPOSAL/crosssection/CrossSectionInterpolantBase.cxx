#include "PROPOSAL/crosssection/CrossSectionInterpolantBase.h"

namespace PROPOSAL {
    Interpolant2DBuilder::Definition CrossSectionInterpolantBase::dNdx_def = {};
    Interpolant1DBuilder::Definition CrossSectionInterpolantBase::dEdx_def = {};
    Interpolant1DBuilder::Definition CrossSectionInterpolantBase::dE2dx_def = {};
} // namespace PROPOSAL
