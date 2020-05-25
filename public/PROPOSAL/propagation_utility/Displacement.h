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

public:
    Displacement(const CrossSectionList&);
    virtual ~Displacement() = default;

    double FunctionToIntegral(const ParticleDef&, const Medium&, double);
    virtual double SolveTrackIntegral( const ParticleDef&, const Medium&, double, double) = 0;
    virtual double UpperLimitTrackIntegral( const ParticleDef&, const Medium&, double, double) = 0;

    size_t GetHash(const ParticleDef&, const Medium&) const;
    double GetLowerLim(const ParticleDef&) const;
};

extern Interpolant1DBuilder::Definition displacement_interpol_def;

template <class T> class DisplacementBuilder : public Displacement {
    unordered_map<size_t, T> integral;
    T BuildTrackIntegral(const ParticleDef&, const Medium&);
public:
    DisplacementBuilder(const CrossSectionList&);
    double SolveTrackIntegral(
        const ParticleDef&, const Medium&, double, double) override;
    double UpperLimitTrackIntegral(
        const ParticleDef&, const Medium&, double, double) override;
};

template <class T>
DisplacementBuilder<T>::DisplacementBuilder(const CrossSectionList& cross)
    : Displacement(cross)
{
}

template <class T>
T DisplacementBuilder<T>::BuildTrackIntegral(
    const ParticleDef& p_def, const Medium& medium)
{
    auto dedx = [this, &p_def, &medium](double energy) {
        return FunctionToIntegral(p_def, medium, energy);
    };
    T integral(dedx, p_def, medium);
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = GetHash(p_def, medium);
        integral.BuildTables("displacement", hash, displacement_interpol_def);
    };
    return integral;
}

template <class T>
double DisplacementBuilder<T>::SolveTrackIntegral(const ParticleDef& p_def,
    const Medium& medium, double upper_lim, double lower_lim)
{
    auto hash = GetHash(p_def, medium);
    auto search = integral.find(hash);
    if (search != integral.end())
        return search->second->Calculate(upper_lim, lower_lim);
    integral[hash] = BuildTrackIntegral(p_def, medium);
    return integral[hash]->Calculate(upper_lim, lower_lim);
}

template <class T>
double DisplacementBuilder<T>::UpperLimitTrackIntegral(const ParticleDef& p_def,
    const Medium& medium, double lower_limit, double sum)
{
    auto hash = GetHash(p_def, medium);
    auto search = integral.find(hash);
    if (search != integral.end())
        return search->second->GetUpperLimit(lower_limit, sum);
    integral[hash] = BuildTrackIntegral(p_def, medium);
    return integral[hash]->GetUpperLimit(lower_limit, sum);
}

} // namespace PROPOSAL
