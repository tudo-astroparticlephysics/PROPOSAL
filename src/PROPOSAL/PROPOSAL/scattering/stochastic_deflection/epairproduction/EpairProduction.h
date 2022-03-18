#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        struct EpairProduction : public Parametrization {
            EpairProduction() = default;
            virtual ~EpairProduction() = default;

            static constexpr InteractionType type = InteractionType::Epair;

            InteractionType GetInteractionType() const noexcept final { return type; };
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL