#pragma once
#include "PROPOSAL/math/InterpolantBuilder.h"

namespace PROPOSAL {
    class Displacement;
}

namespace PROPOSAL {

class Decay {
protected:
    double lifetime;
    double mass;
    std::shared_ptr<Displacement> disp;
    size_t hash;

    double FunctionToIntegral(double energy);

public:
    Decay(std::shared_ptr<Displacement> _disp, double _lifetime, double _mass);
    virtual ~Decay() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double EnergyDecay(double, double, double) = 0;

    auto GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
