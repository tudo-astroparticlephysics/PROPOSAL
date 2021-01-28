#pragma once
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include <cmath>

namespace PROPOSAL {
template <class T> class DecayBuilder : public Decay {
    T decay_integral;

    void build_tables() {};


public:
    DecayBuilder(
        std::shared_ptr<Displacement> _disp, double lifetime, double mass)
        : Decay(_disp, lifetime, mass)
        , decay_integral([this](double E) { return FunctionToIntegral(E); },
              _disp->GetLowerLim(), this->GetHash())
    {
    }

    double EnergyDecay(double energy, double rnd, double density) override
    {
        auto rndd = -std::log(rnd) * density;
        auto rnddMin
            = decay_integral.Calculate(energy, disp->GetLowerLim()) / lifetime;
        if (rndd >= rnddMin)
            return disp->GetLowerLim();
        return decay_integral.GetUpperLimit(energy, rndd * lifetime);
    }
};

inline auto make_decay(std::shared_ptr<Displacement> disp, const ParticleDef& p,
    bool interpolate = true)
{
    auto decay = std::unique_ptr<Decay>();
    if (interpolate)
        decay.reset(
            new DecayBuilder<UtilityInterpolant>(disp, p.lifetime, p.mass));
    else
        decay.reset(
            new DecayBuilder<UtilityIntegral>(disp, p.lifetime, p.mass));
    return decay;
}

template <typename Cross>
inline auto make_decay(
    Cross&& cross, const ParticleDef& p, bool interpolate = true)
{
    auto disp = std::shared_ptr<Displacement>(
        make_displacement(std::forward<Cross>(cross), false));
    return make_decay(disp, p, interpolate);
}
} // namespace PROPOSAL
