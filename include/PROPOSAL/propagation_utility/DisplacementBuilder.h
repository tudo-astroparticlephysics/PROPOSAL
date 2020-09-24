#pragma once
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {
template <typename T>
class DisplacementBuilder : public Displacement {
    T disp_integral;

public:
    T BuildTrackIntegral(crossbase_list_t const& cross)
    {
        if (cross.size() < 1)
            throw std::invalid_argument(
                "at least one crosssection is required.");
        auto disp_func = [this](double energy) {
            return FunctionToIntegral(energy);
        };
        auto low_lim = CrossSectionVector::GetLowerLim(cross);
        T integral(disp_func, low_lim);
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross);
            integral.BuildTables("displacement", hash, interpol_def);
        };
        return integral;
    }

    template <typename Cross>
    DisplacementBuilder(Cross&& cross)
        : Displacement(std::forward<Cross>(cross))
        , disp_integral(BuildTrackIntegral(cross_list))
    {
    }

    inline double SolveTrackIntegral(double upper_lim, double lower_lim) final
    {
        return disp_integral.Calculate(upper_lim, lower_lim);
    }

    inline double UpperLimitTrackIntegral(double lower_lim, double sum) final
    {
        return disp_integral.GetUpperLimit(lower_lim, sum);
    }
};

template <typename T>
std::unique_ptr<Displacement> make_displacement(T&& cross, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<
            DisplacementBuilder<UtilityInterpolant>>(std::forward<T>(cross));
    return PROPOSAL::make_unique<DisplacementBuilder<UtilityIntegral>>(
        std::forward<T>(cross));
}
} // namespace PROPOSAL
