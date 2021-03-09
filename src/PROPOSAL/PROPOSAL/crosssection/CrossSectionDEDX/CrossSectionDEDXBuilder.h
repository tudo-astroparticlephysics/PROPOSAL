#pragma once
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXInterpolant.h"

#include <iostream>

namespace PROPOSAL {
template <typename... Args> auto make_dedx(bool interpolate, Args&&... args)
{
    auto ptr = std::unique_ptr<CrossSectionDEDX>();
    if (interpolate)
        try {
            ptr = std::make_unique<CrossSectionDEDXInterpolant>(
                std::forward<Args>(args)...);
        } catch (exception_axis_builder_dedx_out_of_range const& e) {
            std::cout << e.what() << std::endl;
        }
    else
        ptr = std::make_unique<CrossSectionDEDXIntegral>(
            std::forward<Args>(args)...);
    return ptr;
}
} // namespace PROPOSAL
