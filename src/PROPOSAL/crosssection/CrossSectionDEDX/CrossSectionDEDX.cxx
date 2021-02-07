#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
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

CrossSectionDEDX::CrossSectionDEDX(size_t _hash)
    : hash(_hash)
    , logger(Logging::Get("CrossSection.DEDX"))
{
    logger->info("Creating dEdx.");
}

CrossSectionDEDX::CrossSectionDEDX(crosssection::Parametrization<Medium> const&,
    ParticleDef const&, Medium const&, EnergyCutSettings const&, size_t hash)
    : CrossSectionDEDX(hash)
{
}

CrossSectionDEDX::CrossSectionDEDX(
    crosssection::Parametrization<Component> const&, ParticleDef const&,
    Component const& c, EnergyCutSettings const&, size_t hash)
    : CrossSectionDEDX(detail::generate_dedx_hash(hash, c))
{
}
