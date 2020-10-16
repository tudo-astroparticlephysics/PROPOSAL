#pragma once

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/Parametrization.h"

#include <memory>
#include <map>
#include <vector>

namespace PROPOSAL {
namespace secondaries {

    using TCreateMethod = unique_ptr<Parametrization> (*)(ParticleDef const&, Medium const&);

    class DefaultFactory {
        static std::unique_ptr<std::map<InteractionType, TCreateMethod>> secondaries_map;

    public:
        DefaultFactory() = delete;

        template <typename T>
        static bool Register(InteractionType type)
        {
            if(!secondaries_map)
                secondaries_map = PROPOSAL::make_unique<std::map<InteractionType, TCreateMethod>>();

            auto it = secondaries_map->find(type);
            if (it != secondaries_map->end())
                return false;
            (*secondaries_map)[type] = [](ParticleDef const& p, Medium const& m) {
                return unique_ptr<Parametrization>( PROPOSAL::make_unique<T>(p, m));
            };
            return true;
        }

        static std::unique_ptr<Parametrization> Create(InteractionType type,
                ParticleDef const& p, Medium const& m)
        {
            auto it = secondaries_map->find(type);
            if (it != secondaries_map->end())
                return it->second(p, m);
            std::ostringstream s;
            s << "No secondary builder for this interaction type ("
              << Type_Interaction_Name_Map.find(type)->second << ") available.";
            throw std::logic_error(s.str());
        }
    };

    /* std::map<InteractionType, TCreateMethod> DefaultFactory::secondaries_map = {}; */

} // namespace crosssection
} // namespace PROPOSAL
