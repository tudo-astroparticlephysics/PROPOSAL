#pragma once
#include "PROPOSAL/math/InterpolantBuilder.h"

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


} //namespace PROPOSAL
