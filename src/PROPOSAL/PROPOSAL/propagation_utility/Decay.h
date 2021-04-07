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


public:
    Decay(std::shared_ptr<Displacement> _disp, double _lifetime, double _mass);
    virtual ~Decay() = default;

    virtual double EnergyDecay(double, double, double) = 0;
    double FunctionToIntegral(double energy);

    auto GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
