#pragma once
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

using dndx_map_t = std::unordered_map<std::shared_ptr<const Component>,
    std::unique_ptr<CrossSectionDNDX>>;

template <typename T, typename P1, typename P2, typename P3, typename P4>
auto dndx_builder_cut_dependens(P1 target, P2 param, P3 p_def, P4 cut)
{
    auto dndx = std::unique_ptr<CrossSectionDNDX>();
    if (!cut)
        dndx = std::make_unique<T>(param, p_def, target);
    else
        dndx = std::make_unique<T>(param, p_def, target, *cut);
    return dndx;
}

template <typename... Args> auto dndx_builder(bool interpolate, Args&&... args)
{
    if (interpolate)
        return dndx_builder_cut_dependens<CrossSectionDNDXInterpolant>(
            std::forward<Args>(args)...);
    return dndx_builder_cut_dependens<CrossSectionDNDXIntegral>(
        std::forward<Args>(args)...);
}

template <typename... Args>
dndx_map_t build_cross_section_dndx(
    std::true_type, bool interpol, const Medium& medium, Args&&... args)
{
    auto m = dndx_map_t();
    for (auto& c : medium.GetComponents()) {
        auto comp = std::make_shared<const Components::Component>(c);
        m.emplace(comp, dndx_builder(interpol, c, std::forward<Args>(args)...));
    }
    return m;
}

template <typename... Args>
dndx_map_t build_cross_section_dndx(std::false_type, Args&&... args)
{
    return dndx_map_t { { nullptr },
        { dndx_builder(std::forward<Args>(args)...) } };
}

template <typename Param>
dndx_map_t build_cross_section_dndx(Param param, const ParticleDef& p_def,
    const Medium& medium, std::shared_ptr<const EnergyCutSettings> cut,
    bool interpol)
{
    return build_cross_section_dndx(
        typename Param::base_param_t::component_wise {}, interpol, medium,
        param, p_def, cut);
}
} // namespace PROPOSAL
