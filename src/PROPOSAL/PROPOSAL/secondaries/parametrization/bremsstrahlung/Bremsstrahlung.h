#pragma once

#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Bremsstrahlung : public Parametrization {
        Bremsstrahlung() = default;
        virtual ~Bremsstrahlung() = default;

        static constexpr InteractionType type = InteractionType::Brems;

        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace secondaries
} // namespace PROPOSAL
