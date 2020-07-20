#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"
#include "PROPOSAL/secondaries/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct WeakInteraction : public secondaries::Parametrization {
        WeakInteraction() = default;
        virtual ~WeakInteraction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::WeakInt;
        InteractionType GetInteractionType() const noexcept { return type; };
    };
} // namespace secondaries
} // namespace PROPOSAL
