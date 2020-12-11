#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct EpairProduction : public Parametrization {
        EpairProduction() = default;
        virtual ~EpairProduction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Epair;

        InteractionType GetInteractionType() const noexcept { return type; };

    };
} // namespace secondaries
} // namespace PROPOSAL
