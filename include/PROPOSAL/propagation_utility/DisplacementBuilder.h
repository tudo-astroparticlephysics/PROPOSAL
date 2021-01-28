#pragma once
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

namespace PROPOSAL {
template <typename T> class DisplacementBuilder : public Displacement {
    T disp_integral;

    void build_tables();

public:
    template <typename Cross>
    DisplacementBuilder(Cross&& cross)
        : Displacement(std::forward<Cross>(cross))
        , disp_integral(
              [this](double energy) { return FunctionToIntegral(energy); },
              this->lower_lim, this->hash)
    {
        build_tables();
    }

    inline double SolveTrackIntegral(double lower_lim, double upper_lim) final
    {
        return disp_integral.Calculate(lower_lim, upper_lim);
    }

    inline double UpperLimitTrackIntegral(double lower_lim, double sum) final
    {
        return disp_integral.GetUpperLimit(lower_lim, sum);
    }
};

template <> void DisplacementBuilder<UtilityIntegral>::build_tables();
template <> void DisplacementBuilder<UtilityInterpolant>::build_tables();

template <typename T>
std::unique_ptr<Displacement> make_displacement(T&& cross, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<DisplacementBuilder<UtilityInterpolant>>(
            std::forward<T>(cross));
    return PROPOSAL::make_unique<DisplacementBuilder<UtilityIntegral>>(
        std::forward<T>(cross));
}
} // namespace PROPOSAL
