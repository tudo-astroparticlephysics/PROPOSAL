#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Bremsstrahlung : public secondaries::Parametrization {
        Bremsstrahlung() = default;
        virtual ~Bremsstrahlung() = default;

        static constexpr InteractionType type = PROPOSAL::InteractionType::Brems;

        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace secondaries
} // namespace PROPOSAL
