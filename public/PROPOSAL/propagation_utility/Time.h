#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class Time {
    double mass;

protected:
    template <typename Disp> double FunctionToIntegral(Disp&&, double);

public:
    Time(double mass)
        : mass(mass){};
    virtual ~Time() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double TimeElapsed(double, double, double) = 0;
};

template <typename T, typename Cross> class ExactTimeBuilder : public Time {
    T time_integral;

    T BuildTimeIntegral(Cross&&);

public:
    ExactTimeBuilder(Cross&& cross, const ParticleDef& p_def);

    double TimeElapsed(double, double, double) override;
};

template <typename Disp>
double Time::FunctionToIntegral(Disp&& disp, double energy)
{
    assert(energy > mass);
    auto square_momentum = (energy - mass) * (energy + mass);
    auto particle_momentum = std::sqrt(square_momentum);
    return disp.FunctionToIntegral(energy) * energy
        / (particle_momentum * SPEED);
}

template <typename T, typename Cross>
ExactTimeBuilder<T, Cross>::ExactTimeBuilder(
    Cross&& cross, const ParticleDef& p_def)
    : Time(p_def.mass)
    , time_integral(cross)
{
}

template <typename T, typename Cross>
T ExactTimeBuilder<T, Cross>::BuildTimeIntegral(Cross&& cross)
{
    auto disp = DisplacementBuilder<UtilityIntegral, Cross>(cross);
    auto time_func = [this, &disp](double energy) {
        return FunctionToIntegral(disp, energy);
    };
    T time_integral(time_func, disp.GetLowerLim(cross));
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = disp.GetHash(cross);
        time_integral.BuildTables("time", hash, interpol_def);
    };
    return time_integral;
}

template <typename T, typename Cross>
double ExactTimeBuilder<T, Cross>::TimeElapsed(
    double initial_energy, double final_energy, double)
{
    assert(initial_energy >= final_energy);
    return time_integral.Calculate(initial_energy, final_energy);
}

class ApproximateTimeBuilder : public Time {
public:
    ApproximateTimeBuilder() = default;

    double TimeElapsed(double, double, double distance) override
    {
        assert(distance >= 0);
        return distance / SPEED;
    }
};
}
