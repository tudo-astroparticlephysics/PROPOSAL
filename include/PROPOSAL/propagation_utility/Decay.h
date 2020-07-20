#pragma once
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
#include <cassert>

namespace PROPOSAL {

class Decay {
protected:
    double lifetime;
    double mass;
    double lower_lim;

    template <typename Cross, typename Disp>
    double FunctionToIntegral(Cross&& cross, Disp&& disp, double energy)
    {
        assert(!std::isinf(lifetime));
        assert(energy >= mass);
        double square_momentum = (energy - mass) * (energy + mass);
        double aux = SPEED * std::sqrt(square_momentum) / mass;
        return disp.FunctionToIntegral(cross, energy) / aux;
    }

public:
    Decay(double, double, double);
    virtual ~Decay() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double EnergyDecay(double, double) = 0;
};


} //namespace PROPOSAL
