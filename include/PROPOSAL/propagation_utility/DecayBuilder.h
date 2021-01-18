
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include <math.h>

namespace PROPOSAL {
template <class T> class DecayBuilder : public Decay {
    T decay_integral;

    void build_tables() {};

    /* T BuildDecayIntegral(Cross&& cross) */
    /* { */
    /*     auto disp = std::shared_ptr<Displacement>(make_displacement(cross,
     * false)); */
    /*     auto decay_func = [this, cross, disp](double energy) mutable { */
    /*         return FunctionToIntegral(cross, *disp, energy); */
    /*     }; */
    /*     T decay_integral(decay_func, CrossSectionVector::GetLowerLim(cross));
     */
    /*     if (typeid(T) == typeid(UtilityInterpolant)) { */
    /*         auto hash_digest = (size_t)0; */
    /*         hash_combine(hash_digest, CrossSectionVector::GetHash(cross),
     * interpol_def.GetHash()); */

    /*         decay_integral.BuildTables("decay", hash_digest, interpol_def);
     */
    /*     }; */
    /*     return decay_integral; */
    /* } */

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

auto make_decay(std::shared_ptr<Displacement> disp, const ParticleDef& p_def,
    bool interpolate = true)
{
    auto decay = std::unique_ptr<Decay>();
    if (interpolate)
        decay.reset(new DecayBuilder<UtilityInterpolant>(
            disp, p_def.lifetime, p_def.mass));
    else
        decay.reset(new DecayBuilder<UtilityIntegral>(
            disp, p_def.lifetime, p_def.mass));
    return decay;
}

template <typename Cross>
auto make_decay(
    Cross&& cross, const ParticleDef& p_def, bool interpolate = true)
{
    auto disp = std::shared_ptr<Displacement>(
        make_displacement(std::forward<Cross>(cross), false));
    return make_decay(disp, p_def, interpolate);
}
} // namespace PROPOSAL
