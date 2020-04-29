#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"

using std::shared_ptr;
using std::vector;

namespace PROPOSAL {
class InterpolationDef;
class ParticleDef;

extern InterpolationDef std_interpolation_def;

vector<shared_ptr<CrossSection>> GetStdCrossSections(const Medium&,
    std::shared_ptr<EnergyCutSettings>, const ParticleDef&);

} // namespace PROPOSAL
