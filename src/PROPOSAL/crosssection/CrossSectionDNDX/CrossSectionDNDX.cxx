#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

CrossSectionDNDX::CrossSectionDNDX(lim_func_t _kin_lim,
    std::shared_ptr<const EnergyCutSettings> _cut, size_t _hash)
    : hash(_hash)
    , logger(Logging::Get("CrossSection.DNDX"))
    , kinematic_limits(_kin_lim)
    , cut(_cut)
{
    logger->debug("Creating dNdx.");
}

CrossSectionDNDX::CrossSectionDNDX(param_medium_t const& param, ParticleDef p,
    Medium m, std::shared_ptr<const EnergyCutSettings> cut, size_t hash)
    : CrossSectionDNDX(
        [ptr = std::shared_ptr<param_medium_t>(param.clone()), p, m](
            double E) { return ptr->GetKinematicLimits(p, m, E); },
        cut, hash)
{
}

CrossSectionDNDX::CrossSectionDNDX(param_comp_t const& param, ParticleDef p,
    Component m, std::shared_ptr<const EnergyCutSettings> cut, size_t hash)
    : CrossSectionDNDX(
        [ptr = std::shared_ptr<param_comp_t>(param.clone()), p, m](
            double E) { return ptr->GetKinematicLimits(p, m, E); },
        cut, hash)
{
}

CrossSectionDNDX::IntegrationLimit CrossSectionDNDX::GetIntegrationLimits(
    double energy) const
{
    auto kin_lim = kinematic_limits(energy);
    auto lim = IntegrationLimit { kin_lim.v_min, kin_lim.v_max };
    if (cut)
        lim.min = std::max(cut->GetCut(kin_lim, energy), lim.min);
    return lim;
}
