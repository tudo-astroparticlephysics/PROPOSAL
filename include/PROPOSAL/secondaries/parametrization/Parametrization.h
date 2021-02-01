#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/RegisteredInDefault.h"

#include <vector>

namespace PROPOSAL {
namespace secondaries {
    struct Parametrization {
        Parametrization() = default;
        virtual ~Parametrization() = default;

        virtual size_t RequiredRandomNumbers() const noexcept = 0;
        virtual InteractionType GetInteractionType() const noexcept = 0;
        virtual std::vector<ParticleState> CalculateSecondaries(
            StochasticLoss, const Component&, std::vector<double>& rnd)
            = 0;
    };

    template <typename T>
    using DefaultSecondaries
        = RegisteredInDefault<secondaries::Parametrization, T>;
} // namespace secondaries
} // namespace PROPOSAL
