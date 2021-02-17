#pragma once

#include "PROPOSAL/DefaultFactory.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"
#include "PROPOSAL/scattering/stochastic_deflection/ScatteringFactory.h"
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

    scatter_ptr m_scatter_ptr;
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

protected:
    virtual std::array<double, 2> _scale_deflect(
        std::array<double, 2>& a, InteractionType)
    {
        return a;
    }

    template <typename... Args>
    auto _stochastic_deflect(
        stochastic_deflection::Parametrization const& p, Args... args)
    {
        auto angles = p.CalculateStochasticDeflection(args...);
        return _scale_deflect(angles, p.GetInteractionType());
    }

    virtual std::array<double, 4> _scale_scatter(std::array<double, 4>& a)
    {
        return a;
    }

    template <typename... Args>
    auto _multiple_scatter(
        multiple_scattering::Parametrization& p, Args... args)
    {
        auto angles = p.CalculateRandomAngle(args...);
        return _scale_scatter(angles);
    }

public:
    Scattering() = default;

    /**
     * @brief Storage class of objects related to particle deflection. There are
     * diffentiated between stochastic deflection and multiple scattering
     *
     * @tparam T1 multiple_scattering::Parametrization or nullptr_t
     * @tparam T2 container of stochastic_deflection::Parametrization or
     * nullptr_t
     * @param _m Multiple scattering calculator to take deflections caused by
     * continuous losses into account
     * @param _s list of deflection calculator to take stochastic deflections
     * intor account
     */
    template <typename T1, typename T2>
    Scattering(T1&& _m, T2&& _s)
        : m_scatter_ptr(init_multiple_scatter(std::forward<T1>(_m)))
        , stochastic_deflection(init_deflection(std::forward<T2>(_s)))
    {
    }

    /**
     * @brief random numbers required for a deflection of a specific type.
     */
    size_t StochasticDeflectionRandomNumbers(InteractionType t) const noexcept
    {
        auto it = stochastic_deflection.find(t);
        if (it != stochastic_deflection.end())
            return it->second->RequiredRandomNumbers();
        return 0;
    }

    /**
     * @brief random numbers required for multiple scattering.
     */
    constexpr size_t MultipleScatteringRandomNumbers() noexcept { return 4; }

    /**
     * @brief Calculates deflection angles for specific interaction type. Take a
     * deeper look into stochastic_deflection::Parametrization for a better
     * understanding.
     */
    template <typename... Args>
    std::array<double, 2> CalculateStochasticDeflection(
        InteractionType t, Args... args)
    {
        auto it = stochastic_deflection.find(t);
        if (it != stochastic_deflection.end())
            return _stochastic_deflect(*it->second, args...);
        return std::array<double, 2> { 0., 0. };
    }

    /**
     * @brief Calculate scattering angles in cartesian coordinates. Take a
     * deeper look into multiple_scattering::Parametrization for a better
     * understanding.
     */
    template <typename... Args>
    std::array<double, 4> CalculateMultipleScattering(Args... args)
    {
        if (m_scatter_ptr)
            return _multiple_scatter(*m_scatter_ptr, args...);
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

/**
 * @brief Creates a default deflection of specific type.
 *
 * @tparam Args InteractionType and stochastic_deflection::Parametrization
 * cstr. arguments
 * @param args interaction type and the arguments of the corresponding cstr.
 */
template <typename... Args> inline auto make_default_stochastic_deflection(Args... args)
{
    return DefaultFactory<stochastic_deflection::Parametrization>::Create(
        args...);
}

/**
 * @brief Creates a vector of default deflection for a vector of types.
 *
 * @tparam Args std::vector<InteractionType> and
 * stochastic_deflection::Parametrization cstr. arguments
 * @param types list of interaction types
 * @param args the corresponding cstr. argmuents
 */
template <typename... Args>
inline auto make_default_stochastic_deflection(
    std::vector<InteractionType> const& types, Args... args)
{
    using param_ptr = std::unique_ptr<stochastic_deflection::Parametrization>;
    auto v = std::vector<param_ptr>();
    for (auto t : types)
        v.emplace_back(make_default_stochastic_deflection(t, args...));
    return v;
}
} // namespace PROPOSAL
