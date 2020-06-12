#pragma once

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace secondaries {
    struct Bremsstrahlung {
        Bremsstrahlung() = default;
        virtual ~Bremsstrahlung() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Brems;
    };
} // namespace secondaries
} // namespace PROPOSAL
