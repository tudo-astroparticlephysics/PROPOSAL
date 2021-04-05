#pragma once

#include "PROPOSAL/RegisteredInDefault.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/scattering/DirectionChange.h"

#include <array>
#include <memory>
#include <vector>

using PROPOSAL::Component;

namespace PROPOSAL {
namespace stochastic_deflection {
    struct Parametrization {
        Parametrization() = default;
        virtual ~Parametrization() = default;

        virtual std::unique_ptr<Parametrization> clone() const = 0;

        virtual size_t RequiredRandomNumbers() const noexcept = 0;
        virtual InteractionType GetInteractionType() const noexcept = 0;
        virtual DirectionChangeAngular CalculateStochasticDeflection(
            double initial_energy, double final_energy,
            std::vector<double> const&) const = 0;
    };

    template <typename T>
    using DefaultDeflection
        = RegisteredInDefault<stochastic_deflection::Parametrization, T>;
} // namespace stochastic_deflection
} // namespace PROPOSAL
