#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include <string>
#include <vector>

using std::string;

namespace PROPOSAL {

class Displacement {
    CrossSectionList cross;

public:
    Displacement(const CrossSectionList&);
    virtual ~Displacement() = default;
    virtual double FunctionToIntegral(double energy);
    virtual double SolveTrackIntegral(double, double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;

};

template <class T> class DisplacementBuilder : public Displacement {
    T integral;

public:
    DisplacementBuilder(const CrossSectionList&);
    double SolveTrackIntegral(double, double, double) override;
    double UpperLimitTrackIntegral(double, double) override;

    static Interpolant1DBuilder::Definition displacement_interpol_def;
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
            hash_combine(hash_digest, c->GetHash());
        integral.BuildTables("displacement", hash_digest, displacement_interpol_def);
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

template <class T>
Interpolant1DBuilder::Definition DisplacementBuilder<T>::displacement_interpol_def;
} // namespace PROPOSAL
