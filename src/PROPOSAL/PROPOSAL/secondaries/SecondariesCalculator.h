#pragma once

#include "PROPOSAL/DefaultFactory.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace PROPOSAL {

//!
//! Class to store different secondaries parametrization and produce secondaries
//! for different interactiontypes. Take into account that different particle
//! make use of the same secondary calculator with different parameters.
//!
class SecondariesCalculator {
    using param_ptr = std::unique_ptr<secondaries::Parametrization>;

    //!
    //! Storage of the associated secondary calculator to the corresponding
    //! interaction type.
    //!
    std::unordered_map<InteractionType, param_ptr, InteractionType_hash>
        secondary_generator;

public:
    //!
    //! Empty secondaries calculator. Parametrization has to be added by
    //! addInteraction.
    //!
    SecondariesCalculator() = default;

    //!
    //! Initalize secondary calculator with a bunch of interaction types for a
    //! specific particle and medium.
    //!
    template <typename T>
    SecondariesCalculator(
        T const& interaction_types, ParticleDef const& p, Medium const& m)
    {
        for (auto& t : interaction_types)
            addInteraction(
                DefaultFactory<secondaries::Parametrization>::Create(t, p, m));
    }

    //!
    //! Transfer ownership of a secondary calculator parametrization, to the
    //! secondary calculator.
    //!
    inline void addInteraction(param_ptr&& p)
    {
        secondary_generator[p->GetInteractionType()] = std::move(p);
    }

    //!
    //! Returns number required for calculation of a specific interactiontype.
    //!
    inline size_t RequiredRandomNumbers(InteractionType type) const noexcept
    {
        return secondary_generator.find(type)->second->RequiredRandomNumbers();
    }

    //!
    //! Calculates the secondary particle for a given loss. Initial particle is
    //! treated as a loss and returned as the first particle of the secondaries.
    //!
    std::vector<ParticleState> CalculateSecondaries(
        StochasticLoss, Component const&, std::vector<double>&);
};

//!
//! Produces a SecondariesCalculator.
//!
template <typename TypeList>
std::unique_ptr<SecondariesCalculator> make_secondaries(
    TypeList&& list, ParticleDef const& p, Medium const& m)
{
    return PROPOSAL::make_unique<SecondariesCalculator>(
        std::forward<TypeList>(list), p, m);
}

} // namespace PROPOSAL
