#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"

using std::shared_ptr;
using std::vector;

namespace PROPOSAL {
class InterpolationDef;
class ParticleDef;

extern InterpolationDef std_interpolation_def;

template <typename Param>
shared_ptr<CrossSection> make_crosssection(
    Param&& param, shared_ptr<const EnergyCutSettings> cuts, bool interpolate);

vector<shared_ptr<CrossSection>> GetStdCrossSections(
    std::shared_ptr<EnergyCutSettings>, const ParticleDef&);

} // namespace PROPOSAL
