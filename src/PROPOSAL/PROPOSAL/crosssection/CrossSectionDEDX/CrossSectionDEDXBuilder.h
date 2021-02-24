#pragma once
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXInterpolant.h"

namespace PROPOSAL {

template <typename... Args> auto make_dedx(bool interpolate, Args&&... args)
{
    auto ptr = std::unique_ptr<CrossSectionDEDX>();
    if (interpolate)
        ptr = std::make_unique<CrossSectionDEDXInterpolant>(
            std::forward<Args>(args)...);
    else
        ptr = std::make_unique<CrossSectionDEDXIntegral>(
            std::forward<Args>(args)...);
    return ptr;
}
} // namespace PROPOSAL
