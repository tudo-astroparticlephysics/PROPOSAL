#pragma once
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

using dndx_map_t = std::unordered_map<std::shared_ptr<const Component>,
    std::unique_ptr<CrossSectionDNDX>>;

namespace detail {
    template <typename T, typename P1, typename P2, typename P3, typename P4>
    auto dndx_builder_cut(P1 target, P2 param, P3 p_def, P4 cut)
    {
        if (cut)
            return std::unique_ptr<CrossSectionDNDX>(
                std::make_unique<T>(param, p_def, target, *cut));
        return std::unique_ptr<CrossSectionDNDX>(
            std::make_unique<T>(param, p_def, target));
    }

    template <typename T, typename P1, typename P2, typename P3>
    auto dndx_builder_cut(P1 target, P2 param, P3 p_def, std::nullptr_t)
    {
        return std::unique_ptr<CrossSectionDNDX>(
            std::make_unique<T>(param, p_def, target));
    }

    template <typename... Args>
    auto make_dndx(bool interpolate, Args&&... args)
    {
        if (interpolate)
            return dndx_builder_cut<CrossSectionDNDXInterpolant>(
                std::forward<Args>(args)...);
        return dndx_builder_cut<CrossSectionDNDXIntegral>(
            std::forward<Args>(args)...);
    }
}

template <typename... Args>
inline auto make_dndx(bool interpol, Args&&... args)
{
    return detail::make_dndx(interpol, std::forward<Args>(args)...);
}
} // namespace PROPOSAL
