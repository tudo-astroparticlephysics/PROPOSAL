#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DX.h"
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
    size_t generate_de2dx_hash(size_t hash, Component const& c)
    {
        hash_combine(hash, c.GetHash());
        return hash;
    }
} // namespace detail
} // namespace PROPOSAL

CrossSectionDE2DX::CrossSectionDE2DX(size_t _hash)
    : hash(_hash)
    , logger(Logging::Get("CrossSection.DE2DX"))
{
    logger->info("Creating dEdx.");
}

CrossSectionDE2DX::CrossSectionDE2DX(
    crosssection::Parametrization<Medium> const&, ParticleDef const&,
    Medium const&, EnergyCutSettings const&, size_t hash)
    : CrossSectionDE2DX(hash)
{
}

CrossSectionDE2DX::CrossSectionDE2DX(
    crosssection::Parametrization<Component> const&, ParticleDef const&,
    Component const& c, EnergyCutSettings const&, size_t hash)
    : CrossSectionDE2DX(detail::generate_de2dx_hash(hash, c))
{
}
