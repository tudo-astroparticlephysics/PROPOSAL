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

    template <typename T> inline auto init_deflection(T&& obj)
    {
        auto m = deflect_map_t();
        for (auto&& d : obj)
            m[d->GetInteractionType()] = std::move(d);
        return m;
    }

    template <typename T> inline auto init_deflection(T const& ref)
    {
        auto _copy = std::vector<deflect_ptr>();
        for (auto& _d : ref)
            _copy.emplace_back(_d->clone());
        return init_deflection(std::move(_copy));
    }

    template <typename T> inline auto init_multiple_scatter(T&& obj)
    {
        return std::move(obj);
    }

    template <typename T> inline auto init_multiple_scatter(T const& ref)
    {
        return init_multiple_scatter(ref.clone());
    }

public:
    Scattering() = default;

    template <typename T1, typename T2>
    Scattering(T1&& _m, T2&& _s)
        : multiple_scatter(init_multiple_scatter(std::forward<T1>(_m)))
        , stochastic_deflection(init_deflection(std::forward<T2>(_s)))
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
    std::array<double, 2> CalculateStochasticDeflection(
        InteractionType t, Args... args) const
    {
        auto it = stochastic_deflection.find(static_cast<InteractionType>(t));
        if (it != stochastic_deflection.end())
            return it->second->CalculateStochasticDeflection(args...);
        return std::array<double, 2> { 0., 0. };
    }

    template <typename... Args>
    std::array<double, 4> CalculateMultipleScattering(Args... args) const
    {
        if (multiple_scatter)
            return multiple_scatter->CalculateRandomAngle(args...);
        return std::array<double, 4> { 0, 0, 0, 0 };
    }
};

template <> inline auto Scattering::init_deflection(std::nullptr_t&&)
{
    return Scattering::deflect_map_t();
}

template <> inline auto Scattering::init_multiple_scatter(std::nullptr_t&&)
{
    return nullptr;
}

inline auto make_stochastic_deflection(
    InteractionType t, ParticleDef const& p, Medium const& m)
{
    return DefaultFactory<stochastic_deflection::Parametrization>::Create(
        t, p, m);
}

inline auto make_stochastic_deflection(
    std::vector<InteractionType> const& types, ParticleDef const& p,
    Medium const& m)
{
    using param_ptr = std::unique_ptr<stochastic_deflection::Parametrization>;
    auto v = std::vector<param_ptr>();
    for (auto t : types)
        v.emplace_back(make_stochastic_deflection(t, p, m));
    return v;
}
} // namespace PROPOSAL
