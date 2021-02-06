
#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionInterpolant.h"

namespace PROPOSAL {
template <typename Param, typename P, typename M>
auto make_crosssection(Param&& param, P&& p_def, M&& medium,
    std::shared_ptr<const EnergyCutSettings> cuts, bool interpolate)
{
    auto cross = std::unique_ptr<CrossSectionBase>();
    if (interpolate)
        cross = std::make_unique<CrossSectionInterpolant<Param>>(
            std::forward<Param>(param), std::forward<P>(p_def),
            std::forward<M>(medium), cuts);
    else
        cross = std::make_unique<CrossSectionIntegral<Param>>(
            std::forward<Param>(param), std::forward<P>(p_def),
            std::forward<M>(medium), cuts);
    return cross;
}
}
