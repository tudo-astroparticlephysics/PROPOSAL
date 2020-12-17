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
    using deflect_map_t = std::unordered_map<InteractionType, deflect_ptr,
        InteractionType_hash>;
    using scatter_ptr = std::unique_ptr<multiple_scattering::Parametrization>;

    scatter_ptr multiple_scatter;
    deflect_map_t stochastic_deflection;

    template <typename T> inline deflect_map_t init_deflection(T obj)
    {
        auto m = deflect_map_t();
        for (auto& d : *obj)
            m[d->GetInteractionType()] = d->clone();
        return m;
    }

    template <typename T> inline scatter_ptr init_multiple_scatter(T obj)
    {
        return obj->clone();
    }

public:
    Scattering() = default;

    template <typename T1, typename T2>
    Scattering(T1 _multiple_scatter, T2 _stochastic_deflection)
        : multiple_scatter(init_multiple_scatter(std::move(_multiple_scatter)))
        , stochastic_deflection(init_deflection(std::move(_stochastic_deflection)))
    {
    }

    size_t StochasticDeflectionRandomNumbers(InteractionType t) const noexcept
    {
        auto it = stochastic_deflection.find(t);
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

template <>
inline Scattering::deflect_map_t Scattering::init_deflection(std::nullptr_t)
{
    return Scattering::deflect_map_t();
}

template <>
inline Scattering::scatter_ptr Scattering::init_multiple_scatter(std::nullptr_t)
{
    return nullptr;
}

inline auto make_stochastic_deflection(
    InteractionType t, ParticleDef const& p, Medium const& m)
{
    return DefaultFactory<stochastic_deflection::Parametrization>::Create(
        t, p, m);
}

template <typename T>
inline auto make_stochastic_deflection(
    T types_container, ParticleDef const& p, Medium const& m)
{
    auto v = std::vector<
        std::unique_ptr<stochastic_deflection::Parametrization>>();
    for (auto t : types_container)
        v.emplace_back(make_stochastic_deflection(t, p, m));
    return v;
}
} // namespace PROPOSAL
