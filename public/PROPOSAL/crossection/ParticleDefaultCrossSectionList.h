#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {
class InterpolationDef;
class ParticleDef;

extern InterpolationDef std_interpolation_def;

CrossSectionList GetStdCrossSections(std::shared_ptr<Medium>, std::shared_ptr<EnergyCutSettings>, const ParticleDef&);

} // namespace PROPOSAL
