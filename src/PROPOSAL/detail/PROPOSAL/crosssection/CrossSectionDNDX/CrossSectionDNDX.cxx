#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

CrossSectionDNDX::CrossSectionDNDX(lim_func_t _kin_lim,
    double _lower_energy_lim, std::shared_ptr<const EnergyCutSettings> _cut,
    size_t _hash)
    : hash(_hash)
    , logger(Logging::Get("CrossSection.DNDX"))
    , lower_energy_lim(_lower_energy_lim)
    , kinematic_limits(_kin_lim)
    , cut(_cut)
{
    logger->debug("Creating dNdx.");
}

CrossSectionDNDX::CrossSectionDNDX(param_medium_t const& _param, ParticleDef _p,
    Medium _m, std::shared_ptr<const EnergyCutSettings> _cut, size_t _hash)
    : CrossSectionDNDX(
        [ptr = std::shared_ptr<param_medium_t>(_param.clone()), _p, _m](
            double E) { return ptr->GetKinematicLimits(_p, _m, E); },
        _param.GetLowerEnergyLim(_p), _cut, _hash)
{
}

CrossSectionDNDX::CrossSectionDNDX(param_comp_t const& _param, ParticleDef _p,
    Component _c, std::shared_ptr<const EnergyCutSettings> _cut, size_t _hash)
    : CrossSectionDNDX(
        [ptr = std::shared_ptr<param_comp_t>(_param.clone()), _p, _c](
            double E) { return ptr->GetKinematicLimits(_p, _c, E); },
        _param.GetLowerEnergyLim(_p), _cut, _hash)
{
    hash_combine(this->hash, _c.GetHash());
}

CrossSectionDNDX::IntegrationLimit CrossSectionDNDX::GetIntegrationLimits(
    double energy) const
{
    auto kin_lim = kinematic_limits(energy);
    auto lim = IntegrationLimit();
    lim.min = kin_lim.v_min;
    lim.max = kin_lim.v_max;
    if (cut)
        lim.min = std::max(cut->GetCut(kin_lim, energy), lim.min);
    return lim;
}

double CrossSectionDNDX::GetLowerEnergyLim() const { return lower_energy_lim; }
