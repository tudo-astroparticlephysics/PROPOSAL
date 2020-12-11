#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include <array>
#include <vector>

using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace stochastic_deflection {
    struct Parametrization {
        Parametrization() = default;
        virtual ~Parametrization() = default;

        virtual size_t RequiredRandomNumbers() const noexcept = 0;
        virtual InteractionType GetInteractionType() const noexcept = 0;
        virtual std::array<double, 2> CalculateStochasticDeflection(
            StochasticLoss const&, Component const&, std::vector<double> const&)
            = 0;
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL
