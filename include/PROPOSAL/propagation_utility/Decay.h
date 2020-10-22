#pragma once
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
#include <cassert>

namespace PROPOSAL {

class Decay {
protected:
    double lifetime;
    double mass;
    double lower_lim;

    template <typename Cross>
    double FunctionToIntegral(Cross const&, Displacement & disp, double energy)
    {
        assert(!std::isinf(lifetime));
        assert(energy >= mass);
        double square_momentum = (energy - mass) * (energy + mass);
        double aux = SPEED * std::sqrt(square_momentum) / mass;
        return disp.FunctionToIntegral(energy) / aux;
    }

public:
    Decay(double, double, double);
    virtual ~Decay() = default;

    static std::unique_ptr<Interpolant1DBuilder::Definition> interpol_def;

    virtual double EnergyDecay(double, double, double) = 0;
};


} //namespace PROPOSAL
