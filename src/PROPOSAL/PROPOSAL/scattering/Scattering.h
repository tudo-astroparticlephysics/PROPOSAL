#pragma once

#include "PROPOSAL/DefaultFactory.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"
#include "PROPOSAL/scattering/stochastic_deflection/ScatteringFactory.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace PROPOSAL {
struct ParticleDef;
class Medium;
} // namespace PROPOSAL

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
    virtual DirectionChangeAngular _scale_deflect(
        DirectionChangeAngular& a, InteractionType)
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

    virtual multiple_scattering::ScatteringOffset _scale_scatter(
            multiple_scattering::ScatteringOffset& a)
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
    virtual ~Scattering() = default;

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
    static constexpr size_t MultipleScatteringRandomNumbers() noexcept { return 4; }

    /**
     * @brief Calculates deflection angles for specific interaction type. Take a
     * deeper look into stochastic_deflection::Parametrization for a better
     * understanding.
     */
    template <typename... Args>
    DirectionChangeAngular CalculateStochasticDeflection(
        InteractionType t, Args... args)
    {
        auto it = stochastic_deflection.find(t);
        if (it != stochastic_deflection.end())
            return _stochastic_deflect(*it->second, args...);
        auto new_dir =DirectionChangeAngular();
        new_dir.zenith = 0.;
        new_dir.azimuth = 0.;
        return new_dir;
    }

    /**
     * @brief Calculate scattering angles in cartesian coordinates. Take a
     * deeper look into multiple_scattering::Parametrization for a better
     * understanding.
     */
    template <typename... Args>
    multiple_scattering::ScatteringOffset CalculateMultipleScattering(Args... args)
    {
        if (m_scatter_ptr)
            return _multiple_scatter(*m_scatter_ptr, args...);
        return multiple_scattering::ScatteringOffset();
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

template <typename... Args>
inline auto make_scattering(MultipleScatteringType ms_t,
    std::vector<InteractionType> st_t, ParticleDef const& p, Medium const& m,
    Args... args)
{
    using ms_param = multiple_scattering::Parametrization;
    using st_param = stochastic_deflection::Parametrization;

    auto ms = std::unique_ptr<ms_param>();
    if (ms_t != MultipleScatteringType::NoScattering)
        ms = make_multiple_scattering(ms_t, p, m, std::forward<Args>(args)...);

    auto st = std::vector<std::unique_ptr<st_param>>();
    if (!st_t.empty()) {
        for (auto const& t : st_t)
            st.emplace_back(make_default_stochastic_deflection(t, p, m));
    }

    return std::make_unique<Scattering>(std::move(ms), std::move(st));
}

} // namespace PROPOSAL
