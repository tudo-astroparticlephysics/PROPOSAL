#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        struct PairProduction : public Parametrization {
            PairProduction() = default;
            virtual ~PairProduction() = default;

            static constexpr InteractionType type = InteractionType::Epair;

            InteractionType GetInteractionType() const noexcept final { return type; };
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL