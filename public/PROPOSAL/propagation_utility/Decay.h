#pragma once
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include <math.h>

namespace PROPOSAL {

class Decay {
protected:
    double lifetime;
    double mass;
    double lower_lim;

    template <typename Disp> double FunctionToIntegral(Disp&&, double);

public:
    Decay(double, double, double);
    virtual ~Decay() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double EnergyDecay(double, double) = 0;
};

template <class T, class Cross> class DecayBuilder : public Decay {
    T decay_integral;

    T BuildDecayIntegral(Cross&& cross);

public:
    DecayBuilder<T, Cross>(Cross&&, double, double);

    double EnergyDecay(double, double) override;
};

template <class T, class Cross>
DecayBuilder<T, Cross>::DecayBuilder(
    Cross&& cross, double lifetime, double mass)
    : Decay(lifetime, mass, GetLowerLim(cross))
    , decay_integral(BuildDecayIntegral(cross))
{
}

template <class T, class Cross>
T DecayBuilder<T, Cross>::BuildDecayIntegral(Cross&& cross)
{
    auto disp = DisplacementBuilder<UtilityIntegral, Cross>(cross);
    auto decay_func = [this, &disp](double energy) {
        return FunctionToIntegral(disp, energy);
    };
    T decay_integral(decay_func, disp.GetLowerLim(cross));
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = disp.GetHash(cross);
        decay_integral.BuildTables("decay", hash, interpol_def);
    };
    return decay_integral;
}

template <class T, class Cross>
double DecayBuilder<T, Cross>::EnergyDecay(double initial_energy, double rnd)
{
    auto rndd = -std::log(rnd);
    auto rnddMin
        = decay_integral.Calculate(initial_energy, lower_lim) / lifetime;
    if (rndd >= rnddMin)
        return lower_lim;
    return decay_integral.GetUpperLimit(initial_energy, rndd * lifetime);
}

} //namespace PROPOSAL
