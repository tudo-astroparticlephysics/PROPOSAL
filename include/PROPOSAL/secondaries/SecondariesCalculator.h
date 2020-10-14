#pragma once

#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"
#include "PROPOSAL/secondaries/Parametrization.h"

#include <memory>
#include <unordered_map>
#include <vector>

using std::vector;

namespace PROPOSAL {

class SecondariesCalculator {
    using param_ptr = unique_ptr<secondaries::Parametrization>;

    std::unordered_map<InteractionType, param_ptr, InteractionType_hash> secondary_generator;

public:
    SecondariesCalculator() = default;

    template <typename T>
    SecondariesCalculator(T const& interaction_types, ParticleDef const& p, Medium const& m)
    {
       for (auto& t : interaction_types)
            addInteraction(secondaries::DefaultFactory::Create(t, p, m));
    }

    inline void addInteraction(param_ptr&& p)
    {
        secondary_generator[p->GetInteractionType()] = std::move(p);
    }

    inline size_t RequiredRandomNumbers(InteractionType type) const noexcept
    {
        return secondary_generator.find(type)->second->RequiredRandomNumbers();
    }

    std::vector<Loss::secondary_t> CalculateSecondaries(double, Loss::secondary_t, Component const&, std::vector<double>);
};

template <typename TypeList>
std::unique_ptr<SecondariesCalculator> make_secondaries(TypeList&& list, ParticleDef const& p, Medium const& m)
{
    return PROPOSAL::make_unique<SecondariesCalculator>(std::forward<TypeList>(list), p, m);
}

} // namespace PROPOSAL
