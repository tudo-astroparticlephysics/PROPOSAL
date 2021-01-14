#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        struct NuclearInteraction : public Parametrization {
            NuclearInteraction() = default;
            virtual ~NuclearInteraction() = default;

            static constexpr InteractionType type = InteractionType::Photonuclear;

            InteractionType GetInteractionType() const noexcept final { return type; };
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL