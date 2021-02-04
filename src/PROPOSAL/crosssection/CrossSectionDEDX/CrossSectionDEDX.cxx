#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <iostream>

using namespace PROPOSAL;

CrossSectionDEDX::CrossSectionDEDX(size_t _hash, std::string param_name)
    : hash(_hash)
    , logger(Logging::Get("CrossSection.DEDX"))
{
    hash_combine(hash, param_name);
    logger->info("Building dEdx for parametrization {}", param_name);
}

namespace PROPOSAL {
namespace detail {
    template <typename Target>
    size_t _dEdx_Hash(size_t t,
        crosssection::Parametrization<Target> const& param,
        ParticleDef const& p, Target const& target,
        EnergyCutSettings const& cut)
    {
        size_t hash_diggest;
        hash_combine(hash_diggest, t, param.GetHash(), p.GetHash(),
            target.GetHash(), cut.GetHash());
        return hash_diggest;
    }

    size_t dEdx_Hash(size_t t,
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut)
    {
        return _dEdx_Hash(t, param, p, m, cut);
    }

    size_t dEdx_Hash(size_t t,
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p, Component const& c, EnergyCutSettings const& cut)
    {
        return _dEdx_Hash(t, param, p, c, cut);
    }
} // namespace detail
} // namespace PROPOSAL
