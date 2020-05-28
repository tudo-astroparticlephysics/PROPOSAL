#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

using std::unique_ptr;
using std::vector;

namespace PROPOSAL {

extern InterpolationDef std_interpolation_def;

template <typename Param>
shared_ptr<CrossSection> make_crosssection(
    Param&& param, shared_ptr<const EnergyCutSettings> cuts, bool interpolate);

vector<shared_ptr<CrossSection>> GetStdCrossSections(
    std::shared_ptr<EnergyCutSettings>, const ParticleDef&);

} // namespace PROPOSAL
