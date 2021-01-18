#pragma once
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include <cassert>
#include <cmath>

namespace PROPOSAL {

class Decay {
protected:
    double lifetime;
    double mass;
    std::shared_ptr<Displacement> disp;
    size_t hash;

    auto FunctionToIntegral(double energy)
    {
        assert(!std::isinf(lifetime));
        assert(energy >= mass);
        double square_momentum = (energy - mass) * (energy + mass);
        double aux = SPEED * std::sqrt(square_momentum) / mass;
        return disp->FunctionToIntegral(energy) / aux;
    }

public:
    Decay(std::shared_ptr<Displacement> _disp, double _lifetime, double _mass)
        : lifetime(_lifetime)
        , mass(_mass)
        , disp(_disp)
        , hash(0)
    {
        hash_combine(hash, disp->GetHash(), lifetime, mass);
    }

    virtual ~Decay() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double EnergyDecay(double, double, double) = 0;

    auto GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
