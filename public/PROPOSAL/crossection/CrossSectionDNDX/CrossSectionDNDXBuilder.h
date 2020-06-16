#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

using dndx_map = unordered_map<const Component*, unique_ptr<CrossSectionDNDX>>;

template <typename Param>
dndx_map build_cross_section_dndx(Param param, const ParticleDef& p_def,
    const Medium& medium, shared_ptr<EnergyCutSettings> cut, bool interpol,
    std::true_type)
{
    auto m = dndx_map();
    for (auto const& target : medium.GetComponents()) {
        if (interpol)
            m[&target] = PROPOSAL::make_unique<CrossSectionDNDXInterpolant>(param, p_def, target, cut, InterpolationDef());
        else
            m[&target] = PROPOSAL::make_unique<CrossSectionDNDXIntegral>(param, p_def, target, cut);
    }
    return m;
}

template <typename Param>
dndx_map build_cross_section_dndx(Param param, const ParticleDef& p_def,
    shared_ptr<EnergyCutSettings> cut, const Medium& medium, bool interpol,
    std::false_type)
{
    auto m = dndx_map();
    if (interpol)
        m[nullptr] = PROPOSAL::make_unique<CrossSectionDNDXInterpolant>(param, p_def, medium, cut, InterpolationDef());
    else
        m[nullptr] = PROPOSAL::make_unique<CrossSectionDNDXIntegral>(param, p_def, medium, cut);
    return m;
}

template <typename Param>
dndx_map build_cross_section_dndx(Param param, const ParticleDef& p_def,
    const Medium& medium, shared_ptr<EnergyCutSettings> cut, bool interpol)
{
    return build_cross_section_dndx(param, p_def, medium, cut, interpol,
        typename Param::base_param_t::component_wise{});
}

} // namespace PROPOSAL
