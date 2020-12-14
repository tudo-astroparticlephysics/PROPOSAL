#pragma once

#include "PROPOSAL/RegisteredInDefault.h"
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
            double initial_energy, double final_energy,
            std::vector<double> const&) const = 0;
    };

    template <typename T>
    using DefaultDeflection
        = RegisteredInDefault<stochastic_deflection::Parametrization, T>;
} // namespace stochastic_deflection
} // namespace PROPOSAL
