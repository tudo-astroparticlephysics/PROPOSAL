#pragma once

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/Parametrization.h"

#include <memory>
#include <unordered_map>
#include <vector>

using std::logic_error;
using std::unordered_map;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    using TCreateMethod = unique_ptr<Parametrization> (*)(ParticleDef, Medium);

    class DefaultFactory {
        static unordered_map<InteractionType, TCreateMethod> secondaries_map;

    public:
        DefaultFactory() = delete;

        template <typename T> static bool Register(InteractionType type);
        static unique_ptr<Parametrization> Create(
            const InteractionType& type, ParticleDef p, Medium m);
    };

    template <typename T> bool DefaultFactory::Register(InteractionType type)
    {
        auto it = secondaries_map.find(type);
        if (it == secondaries_map.end()) {
            secondaries_map[type] = [](ParticleDef p, Medium m) {
                return unique_ptr<Parametrization>(
                    PROPOSAL::make_unique<T>(p, m));
            };
            return true;
        }
        std::ostringstream s;
        s << "Not two secondary builder with same interaction type ("
          << Type_Interaction_Name_Map.find(type)->second
          << ") can be added.";
        throw logic_error(s.str());
    }

    template <typename T> class RegisteredInDefault {
    protected:
        static bool s_registered;
        RegisteredInDefault() { s_registered; }
    };

    template <typename T>
    bool RegisteredInDefault<T>::s_registered
        = DefaultFactory::Register<T>(T::GetInteractionType());

} // namespace crosssection
} // namespace PROPOSAL
