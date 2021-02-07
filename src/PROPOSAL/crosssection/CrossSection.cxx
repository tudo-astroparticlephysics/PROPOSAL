#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {
namespace detail {
    template <typename Param>
    size_t _generate_cross_hash(size_t hash, std::string name, unsigned int id,
        Param const& param, ParticleDef const& p, Medium const& m,
        std::shared_ptr<const EnergyCutSettings> cut)
    {
        hash_combine(hash, name, id, param.GetHash(), p.GetHash(), m.GetHash());
        if (cut)
            hash_combine(hash, cut->GetHash());
        return hash;
    }

    size_t generate_cross_hash(size_t hash, std::string name, unsigned int id,
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p, Medium const& m,
        std::shared_ptr<const EnergyCutSettings> cut)
    {
        return _generate_cross_hash(hash, name, id, param, p, m, cut);
    }

    size_t generate_cross_hash(size_t hash, std::string name, unsigned int id,
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p, Medium const& m,
        std::shared_ptr<const EnergyCutSettings> cut)
    {
        return _generate_cross_hash(hash, name, id, param, p, m, cut);
    }
} // namespace detail
} // namespace PROPOSAL
