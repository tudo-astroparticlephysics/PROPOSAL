#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/Time.h"

namespace PROPOSAL {
template <typename T, typename Cross> class ExactTimeBuilder : public Time {
    double mass;
    T time_integral;

    T BuildTimeIntegral(Cross&& cross)
    {
        auto disp = DisplacementBuilder<UtilityIntegral, Cross>(cross);
        auto time_func = [this, &cross, &disp](double energy) {
            return FunctionToIntegral(cross, disp, energy);
        };
        T time_integral(time_func, CrossSectionVector::GetLowerLim(cross));
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross);
            time_integral.BuildTables("time", hash, interpol_def);
        };
        return time_integral;
    }

    template <typename Disp>
    double FunctionToIntegral(Cross&& cross, Disp&& disp, double energy)
    {
        assert(energy > mass);
        auto square_momentum = (energy - mass) * (energy + mass);
        auto particle_momentum = std::sqrt(square_momentum);
        return disp.FunctionToIntegral(cross, energy) * energy
            / (particle_momentum * SPEED);
    }

public:
    ExactTimeBuilder(Cross&& cross, const ParticleDef& p_def)
        : mass(p_def.mass)
        , time_integral(BuildTimeIntegral(cross))
    {
    }

    double TimeElapsed(
        double initial_energy, double final_energy, double) override
    {
        assert(initial_energy >= final_energy);
        return time_integral.Calculate(initial_energy, final_energy);
    }
};


template <typename T>
std::unique_ptr<Time> make_time(
    T&& cross, const ParticleDef& p_def, bool interpolate = true)
{
    if (interpolate)
        return PROPOSAL::make_unique<ExactTimeBuilder<UtilityInterpolant, T>>(
            std::forward<T>(cross), p_def);
    return PROPOSAL::make_unique<ExactTimeBuilder<UtilityIntegral, T>>(
            std::forward<T>(cross), p_def);
}

struct ApproximateTimeBuilder : public Time {
    ApproximateTimeBuilder() = default;

    double TimeElapsed(double, double, double distance) override
    {
        assert(distance >= 0);
        return distance / SPEED;
    }
};
} // namespace PROPOSAL
