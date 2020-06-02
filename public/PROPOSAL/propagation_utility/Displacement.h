#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <string>
#include <vector>

using std::string;

namespace PROPOSAL {

struct Displacement {
    Displacement() = default;
    virtual ~Displacement() = default;

    template <typename Cross> double FunctionToIntegral(Cross&&, double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;

protected:
    template <typename Cross> size_t GetHash(Cross&&) const;
    template <typename Cross> double GetLowerLim(Cross&&) const;
};

extern Interpolant1DBuilder::Definition displacement_interpol_def;

template <class T, class Cross>
class DisplacementBuilder : public Displacement {
    Cross crosssection_list;
    T disp_integral;

protected:
    T BuildTrackIntegral();

public:
    DisplacementBuilder(Cross&&);
    double SolveTrackIntegral(double, double) override;
    double UpperLimitTrackIntegral(double, double) override;
};

template <class T, class Cross>
DisplacementBuilder<T, Cross>::DisplacementBuilder(Cross&& cross)
    : crosssection_list(cross)
    , disp_integral(BuildTrackIntegral())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");
}

template <class T, class Cross>
T DisplacementBuilder<T, Cross>::BuildTrackIntegral()
{
    auto dedx = [this](double energy) {
        return FunctionToIntegral(crosssection_list, energy);
    };
    T integral(dedx, GetLowerLim(crosssection_list));
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = GetHash(crosssection_list);
        integral.BuildTables("displacement", hash, displacement_interpol_def);
    };
    return integral;
}

template <class T, class Cross>
double DisplacementBuilder<T, Cross>::SolveTrackIntegral(
    double upper_lim, double lower_lim)
{
    return disp_integral->Calculate(upper_lim, lower_lim);
}

template <class T, class Cross>
double DisplacementBuilder<T, Cross>::UpperLimitTrackIntegral(
    double lower_limit, double sum)
{
    return disp_integral->GetUpperLimit(lower_limit, sum);
}

} // namespace PROPOSAL
