#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <string>
#include <vector>
#include <utility>

using std::string;

namespace PROPOSAL {

struct Displacement {
    static Interpolant1DBuilder::Definition interpol_def;

    Displacement() = default;
    virtual ~Displacement() = default;

    template <typename Cross> double FunctionToIntegral(Cross&&, double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;
};

template <typename T, typename Cross>
class DisplacementBuilder : public Displacement {
    T disp_integral;

    T BuildTrackIntegral(const Cross& cross);

public:
    DisplacementBuilder(Cross&&);
    double SolveTrackIntegral(double, double) override;
    double UpperLimitTrackIntegral(double, double) override;
};

template <typename Cross>
double Displacement::FunctionToIntegral(Cross&& cross, double energy)
{
    auto result = 0.0;
    for (auto& cr : cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}

struct CrossSectionVector {
    template <typename CrossVec> static double GetLowerLim(CrossVec&& cross_vec)
    {
        using val_t = typename std::decay<CrossVec>::type::value_type;
        auto result = std::max_element(
            cross_vec.begin(), cross_vec.end(), [](val_t a, val_t b) {
                return a->GetLowerEnergyLim() < b->GetLowerEnergyLim();
            });
        return (*result)->GetLowerEnergyLim();
    }

    template <typename CrossVec> static size_t GetHash(CrossVec&& cross)
    {
        auto hash_digest = size_t{ 0 };
        for (auto& c : cross)
            hash_combine(hash_digest, c->GetHash());
        return hash_digest;
    }
};

template <typename T, typename Cross>
DisplacementBuilder<T, Cross>::DisplacementBuilder(Cross&& cross)
    : Displacement()
    , disp_integral(BuildTrackIntegral(std::forward<Cross>(cross)))
{
}

template <typename T, typename Cross>
T DisplacementBuilder<T, Cross>::BuildTrackIntegral(const Cross& cross)
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");
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

template <typename T, typename Cross>
double DisplacementBuilder<T, Cross>::SolveTrackIntegral(
    double upper_lim, double lower_lim)
{
    return disp_integral.Calculate(upper_lim, lower_lim);
}

template <typename T, typename Cross>
double DisplacementBuilder<T, Cross>::UpperLimitTrackIntegral(
    double lower_limit, double sum)
{
    return disp_integral.GetUpperLimit(lower_limit, sum);
}

} // namespace PROPOSAL
