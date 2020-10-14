#pragma once

#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/Constants.h"
namespace PROPOSAL {
template <typename T> class ExactTimeBuilder : public Time {
    double mass;
    T time_integral;

    template <typename Cross>
    T BuildTimeIntegral(Cross const& cross)
    {
        auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
        auto time_func = [this, cross, disp](double energy)  {
            return FunctionToIntegral(cross, *disp, energy);
        };
        T time_integral(time_func, CrossSectionVector::GetLowerLim(cross));
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross);
            time_integral.BuildTables("time", hash, interpol_def);
        };
        return time_integral;
    }

    template <typename Cross>
    double FunctionToIntegral(Cross const&, Displacement& disp, double energy)
    {
        assert(energy > mass);
        auto square_momentum = (energy - mass) * (energy + mass);
        auto particle_momentum = std::sqrt(square_momentum);
        return disp.FunctionToIntegral(energy) * energy
            / (particle_momentum * SPEED);
    }

public:
    template <typename Cross>
    ExactTimeBuilder(Cross&& cross, const ParticleDef& p_def)
        : mass(p_def.mass)
        , time_integral(BuildTimeIntegral(std::forward<Cross>(cross)))
    {
    }

    double TimeElapsed(
        double initial_energy, double final_energy, double distance, double density) override
    {
        (void)distance;
        assert(initial_energy >= final_energy);
        return time_integral.Calculate(initial_energy, final_energy) / density;
    }
};


template <typename T>
std::unique_ptr<Time> make_time(
    T&& cross, const ParticleDef& p_def, bool interpolate = true)
{
    if (interpolate)
        return PROPOSAL::make_unique<ExactTimeBuilder<UtilityInterpolant>>(
            std::forward<T>(cross), p_def);
    return PROPOSAL::make_unique<ExactTimeBuilder<UtilityIntegral>>(
            std::forward<T>(cross), p_def);
}

struct ApproximateTimeBuilder : public Time {
    ApproximateTimeBuilder() = default;

    double TimeElapsed(double, double, double distance, double density) override
    {
        (void) density;
        assert(distance >= 0);
        return distance / SPEED;
    }
};
} // namespace PROPOSAL
