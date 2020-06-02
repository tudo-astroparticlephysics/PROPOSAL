#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

using std::make_shared;
using std::unique_ptr;
using std::vector;

namespace PROPOSAL {

extern InterpolationDef std_interpolation_def;

template <typename Param, typename P, typename M>
shared_ptr<CrossSection<P, M>> make_crosssection(
    Param&& param, P&& p_def, M&& medium, shared_ptr<const EnergyCutSettings> cuts, bool interpolate){
    if (interpolate)
        return make_shared<CrossSectionInterpolant<Param, P, M>>(
            param, p_def, medium, cuts, InterpolationDef());
    return make_shared<CrossSectionIntegral<Param, P, M>>(
        param, p_def, medium, cuts);
}

/* template <typename P, typename M> */
/* vector<shared_ptr<CrossSection<P, M>>> GetStdCrossSections(P&&, M&&, std::shared_ptr<EnergyCutSettings>); */

} // namespace PROPOSAL
