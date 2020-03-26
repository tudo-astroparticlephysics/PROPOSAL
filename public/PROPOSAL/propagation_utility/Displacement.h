#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <string>
#include <vector>

using std::string;

namespace PROPOSAL {
class Displacement {
    CrossSectionList cross;

public:
    Displacement(CrossSectionList);
    virtual ~Displacement() = default;
    virtual double FunctionToIntegral(double energy);
    virtual double SolveTrackIntegral(double, double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;

    string name{ "displacement" };
};

template <class T> class DisplacementBuilder : public Displacement {
    T integral;

public:
    DisplacementBuilder(const CrossSectionList&);
    double SolveTrackIntegral(double, double, double) override;
    double UpperLimitTrackIntegral(double, double) override;
};

template <class T>
DisplacementBuilder<T>::DisplacementBuilder(const CrossSectionList& cross)
    : Displacement(cross)
    , integral(std::bind(
          &Displacement::FunctionToIntegral, this, std::placeholders::_1))
{
    if (typeid(T) == typeid(UtilityInterpolant)) {
        size_t hash_digest{ 0 };
        for (const auto& c : cross)
            hash_combine(hash_digest, c->GetParametrization().GetHash(),
                c->GetParametrization().GetMultiplier());
        integral.BuildTables(name, hash_digest);
    }
}

template <class T>
double DisplacementBuilder<T>::SolveTrackIntegral(
    double upper_limit, double lower_limit, double sum)
{
    return integral.Calculate(upper_limit, lower_limit, sum);
}

template <class T>
double DisplacementBuilder<T>::UpperLimitTrackIntegral(
    double lower_limit, double sum)
{
    return integral.GetUpperLimit(lower_limit, sum);
}
} // namespace PROPOSAL
