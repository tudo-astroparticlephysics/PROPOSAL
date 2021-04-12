#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/particle/Particle.h"

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

    double calculate_lower_energy_lim(
        std::vector<std::tuple<double, std::unique_ptr<CrossSectionDEDX>>>*
            dedx_ptr)
    {
        if (!dedx_ptr)
            return 0.;
        auto min = std::numeric_limits<double>::infinity();
        for (auto& it : *dedx_ptr) {
            auto n = std::get<1>(it)->GetLowerEnergyLim();
            min = (min < n) ? min : n;
        }
        return min;
    }

    std::shared_ptr<spdlog::logger> init_logger(std::string const& param_name,
        size_t id, ParticleDef const& p, Medium const& m,
        std::shared_ptr<const EnergyCutSettings> cut)
    {
        auto logger = Logging::Get("Crosssection");
        auto interaction_type =static_cast<InteractionType>(id);
        auto interaction_name = Type_Interaction_Name_Map.at(interaction_type);
        logger->info("Building {} cross section for interaction type {}.", param_name,
            interaction_name);
        logger->debug("-> with particle {}", p.name);
        logger->debug("-> with target {}", m.GetName());
        if (cut)
            logger->debug("-> with e cut {} and v_cut {}", cut->GetEcut(),
                cut->GetVcut());
        return logger;
    };
} // namespace detail
} // namespace PROPOSAL
