#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
    namespace secondaries {
        struct Photoeffect : public secondaries::Parametrization {
            Photoeffect() = default;
            virtual ~Photoeffect() = default;

            static constexpr InteractionType type = PROPOSAL::InteractionType::Photoeffect;
            InteractionType GetInteractionType() const noexcept { return type; };

        };
    } // namespace secondaries
} // namespace PROPOSAL
