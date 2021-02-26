#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

namespace PROPOSAL {
namespace detail {
    size_t generate_dedx_hash(size_t hash, Component const& c)
    {
        hash_combine(hash, c.GetHash());
        return hash;
    }
} // namespace detail
} // namespace PROPOSAL

CrossSectionDEDX::CrossSectionDEDX(double _lower_energy_lim, size_t _hash)
    : hash(_hash)
    , logger(Logging::Get("CrossSection.DEDX"))
    , lower_energy_lim(_lower_energy_lim)
{
    logger->debug("Creating dEdx.");
}

CrossSectionDEDX::CrossSectionDEDX(
    crosssection::Parametrization<Medium> const& param, ParticleDef const& p,
    Medium const&, EnergyCutSettings const&, size_t hash)
    : CrossSectionDEDX(param.GetLowerEnergyLim(p), hash)
{
}

CrossSectionDEDX::CrossSectionDEDX(
    crosssection::Parametrization<Component> const& param, ParticleDef const& p,
    Component const& c, EnergyCutSettings const&, size_t hash)
    : CrossSectionDEDX(
        param.GetLowerEnergyLim(p), detail::generate_dedx_hash(hash, c))
{
}
