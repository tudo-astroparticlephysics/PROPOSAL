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
            auto hash_digest = (size_t)0;
            hash_combine(hash_digest, CrossSectionVector::GetHash(cross), interpol_def.GetHash());

            time_integral.BuildTables("time", hash_digest, interpol_def);
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
        double initial_energy, double final_energy, double grammage, double local_density) override
    {
        (void)grammage;
        assert(initial_energy >= final_energy);
        return time_integral.Calculate(initial_energy, final_energy) / local_density;
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

    double TimeElapsed(double, double, double grammage, double local_density) override
    {
        assert(grammage >= 0);
        return grammage / (local_density * SPEED);
    }
};
} // namespace PROPOSAL
