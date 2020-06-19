#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"
#include "PROPOSAL/secondaries/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct WeakInteraction : public secondaries::Parametrization {
        WeakInteraction() = default;
        virtual ~WeakInteraction() = default;

        static InteractionType GetInteractionType()
        {
            return PROPOSAL::InteractionType::WeakInt;
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
