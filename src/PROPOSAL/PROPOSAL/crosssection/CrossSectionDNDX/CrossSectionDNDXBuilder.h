#pragma once
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"

namespace PROPOSAL {

using dndx_map_t = std::unordered_map<std::shared_ptr<const Component>,
    std::unique_ptr<CrossSectionDNDX>>;

template <typename... Args> auto make_dndx(bool interpolate, Args&&... args)
{
    auto dndx = std::unique_ptr<CrossSectionDNDX>();
    if (interpolate)
        dndx = std::make_unique<CrossSectionDNDXInterpolant>(
            std::forward<Args>(args)...);
    else
        dndx = std::make_unique<CrossSectionDNDXIntegral>(
            std::forward<Args>(args)...);
    return dndx;
}
} // namespace PROPOSAL
