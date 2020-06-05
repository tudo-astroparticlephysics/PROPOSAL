
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include <math.h>

namespace PROPOSAL {
template <class T, class Cross> class DecayBuilder : public Decay {
    T decay_integral;

    T BuildDecayIntegral(Cross&& cross)
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

public:
    DecayBuilder<T, Cross>(Cross&& cross, double lifetime, double, double mass)
        : Decay(lifetime, mass, GetLowerLim(cross))
        , decay_integral(BuildDecayIntegral(cross))
    {
    }

    double EnergyDecay(double energy, double rnd) override
    {
        auto rndd = -std::log(rnd);
        auto rnddMin = decay_integral.Calculate(energy, lower_lim) / lifetime;
        if (rndd >= rnddMin)
            return lower_lim;
        return decay_integral.GetUpperLimit(energy, rndd * lifetime);
    }
};

template <typename T>
std::unique_ptr<Decay> make_decay(
    T&& cross, double lifetime, double mass, bool interpolate = true)
{
    if (interpolate)
        return PROPOSAL::make_unique<DecayBuilder<UtilityInterpolant, T>>(
            std::forward<T>(cross), lifetime, mass);
    return PROPOSAL::make_unique<DecayBuilder<UtilityIntegral, T>>(
            std::forward<T>(cross), lifetime, mass);
}
} // namespace PROPOSAL
