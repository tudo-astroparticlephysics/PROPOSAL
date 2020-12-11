#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace stochastic_deflection {
    struct Bremsstrahlung : public stochastic_deflection::Parametrization {
        Bremsstrahlung() = default;
        virtual ~Bremsstrahlung() = default;

        static constexpr InteractionType type = PROPOSAL::InteractionType::Brems;

        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL
