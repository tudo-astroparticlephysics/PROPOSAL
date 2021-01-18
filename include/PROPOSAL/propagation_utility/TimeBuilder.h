#pragma once

#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/Time.h"
namespace PROPOSAL {
template <typename T> class ExactTimeBuilder : public Time {
    T time_integral;

public:
    ExactTimeBuilder(std::shared_ptr<Displacement> _disp, double _mass)
        : Time(_disp, _mass)
        , time_integral([this](double E) { return FunctionToIntegral(E); },
              _disp->GetLowerLim(), this->GetHash())
    {
    }

    double TimeElapsed(double initial_energy, double final_energy,
        double grammage, double local_density) override
    {
        (void)grammage;
        assert(initial_energy >= final_energy);
        return time_integral.Calculate(initial_energy, final_energy)
            / local_density;
    }
};

inline auto make_time(std::shared_ptr<Displacement> disp, const ParticleDef& p_def,
    bool interpolate = true)
{
    auto time = std::unique_ptr<Time>();
    if (interpolate)
        time.reset(new ExactTimeBuilder<UtilityInterpolant>(disp, p_def.mass));
    else
        time.reset(new ExactTimeBuilder<UtilityIntegral>(disp, p_def.mass));
    return time;
}

template <typename T>
auto make_time(T&& cross, const ParticleDef& p_def, bool interpolate = true)
{
    auto disp = std::shared_ptr<Displacement>(
        make_displacement(std::forward<T>(cross), false));
    return make_time(disp, p_def, interpolate);
}

struct ApproximateTimeBuilder : public Time {
    ApproximateTimeBuilder() = default;

    template <typename... Args> ApproximateTimeBuilder(Args...) { }

    ~ApproximateTimeBuilder() = default;

    double TimeElapsed(
        double, double, double grammage, double local_density) override
    {
        assert(grammage >= 0);
        return grammage / (local_density * SPEED);
    }
};
} // namespace PROPOSAL
