#pragma once

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace secondaries {
    struct WeakInteraction {
        WeakInteraction() = default;
        virtual ~WeakInteraction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::WeakInt;
    };
} // namespace secondaries
} // namespace PROPOSAL
