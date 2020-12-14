#pragma once

#include "PROPOSAL/DefaultFactory.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"
#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace PROPOSAL {

class Scattering {
    using deflect_ptr = std::unique_ptr<stochastic_deflection::Parametrization>;

    std::unordered_map<InteractionType, deflect_ptr, InteractionType_hash>
        stochastic_deflection;
    std::unique_ptr<multiple_scattering::Parametrization> multiple_scatter;

public:
    Scattering() = default;

    Scattering(
        std::unique_ptr<multiple_scattering::Parametrization> _multiple_scatter,
        std::unique_ptr<std::vector<deflect_ptr>> _stochastic_deflection)
        : multiple_scatter(std::move(_multiple_scatter))
    {
        if (_stochastic_deflection) {
            for (auto& d : *_stochastic_deflection) {
                stochastic_deflection[d->GetInteractionType()] = std::move(d);
            }
        }
    }

    size_t StochasticDeflectionRandomNumbers(
        InteractionType type) const noexcept
    {
        auto it = stochastic_deflection.find(type);
        if (it != stochastic_deflection.end())
            return it->second->RequiredRandomNumbers();
        return 0;
    }

    constexpr size_t MultipleScatteringRandomNumbers() noexcept { return 4; }

    template <typename... Args>
    auto CalculateStoachsticDeflection(InteractionType type, Args... args) const
    {
        auto it
            = stochastic_deflection.find(static_cast<InteractionType>(type));
        if (it != stochastic_deflection.end())
            return it->second->CalculateStochasticDeflection(args...);
        return std::array<double, 2> { 0., 0. };
    }

    template <typename... Args>
    auto CalculateMultipleScattering(Args... args) const
    {
        return multiple_scatter->Scatter(args...);
    }
};

std::unique_ptr<std::vector<std::unique_ptr<stochastic_deflection::Parametrization>>>
make_stochastic_deflection(std::vector<InteractionType> const& types,
    ParticleDef const& p, Medium const& m);
} // namespace PROPOSAL
