#pragma once
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {
template <typename T, typename Cross>
class DisplacementBuilder : public Displacement {
    T disp_integral;

    T BuildTrackIntegral(const Cross& cross)
    {
        if (cross.size() < 1)
            throw std::invalid_argument(
                "at least one crosssection is required.");
        auto disp_func = [this, &cross](double energy) {
            return FunctionToIntegral(cross, energy);
        };
        auto low_lim = CrossSectionVector::GetLowerLim(cross);
        T integral(disp_func, low_lim);
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross);
            integral.BuildTables("displacement", hash, interpol_def);
        };
        return integral;
    }

public:
    DisplacementBuilder(Cross&& cross)
        : Displacement()
        , disp_integral(BuildTrackIntegral(std::forward<Cross>(cross)))
    {
    }

    double SolveTrackIntegral(double upper_lim, double lower_lim) override
    {
        return disp_integral.Calculate(upper_lim, lower_lim);
    }

    double UpperLimitTrackIntegral(double lower_lim, double sum) override
    {
        return disp_integral.GetUpperLimit(lower_lim, sum);
    }
};

template <typename T>
std::unique_ptr<Displacement> make_displacement(T&& cross, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<
            DisplacementBuilder<UtilityInterpolant, T>>(std::forward<T>(cross));
    return PROPOSAL::make_unique<DisplacementBuilder<UtilityIntegral, T>>(
        std::forward<T>(cross));
}
} // namespace PROPOSAL
