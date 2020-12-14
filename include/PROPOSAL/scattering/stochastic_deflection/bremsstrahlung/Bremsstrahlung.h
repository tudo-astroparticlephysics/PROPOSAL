#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"

namespace PROPOSAL {
namespace stochastic_deflection {
    struct Bremsstrahlung : public Parametrization {
        Bremsstrahlung() = default;
        virtual ~Bremsstrahlung() = default;

        static constexpr InteractionType type = InteractionType::Brems;

        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL
