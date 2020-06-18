#pragma once

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace secondaries {
    struct Photonuclear {
        Photonuclear() = default;
        virtual ~Photonuclear() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photonuclear;
    };
} // namespace secondaries
} // namespace PROPOSAL
