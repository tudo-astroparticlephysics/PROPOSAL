#pragma once
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

using dndx_map_t = unordered_map<std::shared_ptr<const Component>, unique_ptr<CrossSectionDNDX>>;

template <typename Param>
dndx_map_t build_cross_section_dndx(Param param, const ParticleDef& p_def,
    const Medium& medium, shared_ptr<const EnergyCutSettings> cut, bool interpol,
    std::true_type)
{
    auto m = dndx_map_t();
    for (auto target : medium.GetComponents()) {
        if (interpol)
            m.emplace(
                    std::make_shared<const Components::Component>(target),
                    PROPOSAL::make_unique<CrossSectionDNDXInterpolant>(param, p_def, target, cut, InterpolationDef())
            );
        else
            m.emplace(
                    std::make_shared<const Components::Component>(target),
                    PROPOSAL::make_unique<CrossSectionDNDXIntegral>(param, p_def, target, cut)
            );
    }
    return m;
}

template <typename Param>
dndx_map_t build_cross_section_dndx(Param param, const ParticleDef& p_def,
    const Medium& medium, shared_ptr<const EnergyCutSettings> cut,  bool interpol,
    std::false_type)
{
    auto m = dndx_map_t();
    if (interpol)
        m.emplace(nullptr, PROPOSAL::make_unique<CrossSectionDNDXInterpolant>(param, p_def, medium, cut, InterpolationDef()));
    else
        m.emplace(nullptr, PROPOSAL::make_unique<CrossSectionDNDXIntegral>(param, p_def, medium, cut));
    return m;
}

template <typename Param>
dndx_map_t build_cross_section_dndx(Param param, const ParticleDef& p_def,
    const Medium& medium, shared_ptr<const EnergyCutSettings> cut, bool interpol)
{
    return build_cross_section_dndx(param, p_def, medium, cut, interpol,
        typename Param::base_param_t::component_wise{});
}

} // namespace PROPOSAL
