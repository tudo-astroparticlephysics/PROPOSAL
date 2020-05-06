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

class Displacement {
    CrossSectionList cross;

protected:
    double lower_lim;

public:
    Displacement(const CrossSectionList&);
    virtual ~Displacement() = default;

    double FunctionToIntegral(double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;
};

extern Interpolant1DBuilder::Definition displacement_interpol_def;

template <class T> class DisplacementBuilder : public Displacement {
    T integral;

public:
    DisplacementBuilder(const CrossSectionList&);
    double SolveTrackIntegral(double, double) override;
    double UpperLimitTrackIntegral(double, double) override;
};

template <class T>
DisplacementBuilder<T>::DisplacementBuilder(const CrossSectionList& cross)
    : Displacement(cross)
    , integral(std::bind(&Displacement::FunctionToIntegral, this,
                   std::placeholders::_1),
          lower_lim)
{

    if (typeid(T) == typeid(UtilityInterpolant)) {
        size_t hash_digest{ 0 };
        for (const auto& c : cross)
            hash_combine(hash_digest, c->GetHash());
        displacement_interpol_def.function1d = [this](double energy) {
            return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate(
                energy, lower_lim);
        };
        integral.BuildTables(
            "displacement", hash_digest, displacement_interpol_def);
    }
}

template <class T>
double DisplacementBuilder<T>::SolveTrackIntegral(
    double upper_limit, double lower_limit)
{
    return integral.Calculate(upper_limit, lower_limit);
}

template <class T>
double DisplacementBuilder<T>::UpperLimitTrackIntegral(
    double lower_limit, double sum)
{
    return integral.GetUpperLimit(lower_limit, sum);
}

} // namespace PROPOSAL
