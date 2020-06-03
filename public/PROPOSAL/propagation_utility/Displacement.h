#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <string>
#include <vector>

using std::string;

namespace PROPOSAL {

struct Displacement {
    static Interpolant1DBuilder::Definition interpol_def;

    Displacement() = default;
    virtual ~Displacement() = default;

    template <typename Cross> double FunctionToIntegral(Cross&&, double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;

protected:
    template <typename Cross> size_t GetHash(Cross&&) const;
};

template <class T, class Cross>
class DisplacementBuilder : public Displacement {
    T disp_integral;

    T BuildTrackIntegral(Cross&&);

public:
    DisplacementBuilder(Cross&&);
    double SolveTrackIntegral(double, double) override;
    double UpperLimitTrackIntegral(double, double) override;
};

template <typename Cross> double GetLowerLim(Cross&& cross)
{
    auto result
        = std::max_element(cross.begin(), cross.end(), [](Cross a, Cross b) {
              return b->GetLowerEnergyLim() < b->GetLowerEnergyLim();
          });
    return *result->GetLowerEnergyLim();
}

template <typename Cross> size_t GetHash(Cross&& cross)
{
    auto hash_digest = size_t{ 0 };
    for (const auto& c : cross)
        hash_combine(hash_digest, c->GetHash());
    return hash_digest;
}

template <class T, class Cross>
DisplacementBuilder<T, Cross>::DisplacementBuilder(Cross&& cross)
    : disp_integral(BuildTrackIntegral())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");
}

template <class T, class Cross>
T DisplacementBuilder<T, Cross>::BuildTrackIntegral(Cross&& cross)
{
    auto disp_func = [this, &cross](double energy) {
        return FunctionToIntegral(cross, energy);
    };
    T integral(disp_func, GetLowerLim(cross));
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = GetHash(cross);
        integral.BuildTables("displacement", hash, interpol_def);
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
