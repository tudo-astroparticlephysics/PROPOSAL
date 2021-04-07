#pragma once
#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXInterpolant.h"

namespace PROPOSAL {

template <typename... Args> auto make_de2dx(bool interpolate, Args&&... args)
{
    auto ptr = std::unique_ptr<CrossSectionDE2DX>();
    if (interpolate)
        ptr = std::make_unique<CrossSectionDE2DXInterpolant>(
            std::forward<Args>(args)...);
    else
        ptr = std::make_unique<CrossSectionDE2DXIntegral>(
            std::forward<Args>(args)...);
    return ptr;
}
} // namespace PROPOSAL
