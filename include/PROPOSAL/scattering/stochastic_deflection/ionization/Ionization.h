#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        struct Ionization : public Parametrization {
            Ionization() = default;
            virtual ~Ionization() = default;

            static constexpr InteractionType type = InteractionType::Ioniz;

            InteractionType GetInteractionType() const noexcept final { return type; };
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL
