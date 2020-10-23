
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include <math.h>

namespace PROPOSAL {
template <class T, class Cross> class DecayBuilder : public Decay {
    T decay_integral;

    T BuildDecayIntegral(Cross&& cross)
    {
        auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
        auto decay_func = [this, cross, disp](double energy) mutable {
            return FunctionToIntegral(cross, *disp, energy);
        };
        T decay_integral(decay_func, CrossSectionVector::GetLowerLim(cross));
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash_digest = (size_t)0;
            hash_combine(hash_digest, CrossSectionVector::GetHash(cross), interpol_def.GetHash());

            decay_integral.BuildTables("decay", hash_digest, interpol_def);
        };
        return decay_integral;
    }

public:
    DecayBuilder<T, Cross>(Cross&& cross, double lifetime, double mass)
        : Decay(lifetime, mass, CrossSectionVector::GetLowerLim(cross))
        , decay_integral(BuildDecayIntegral(cross))
    {
    }

    DecayBuilder<T, Cross>(Cross&& cross, const ParticleDef& p_def)
        : DecayBuilder(cross, p_def.lifetime, p_def.mass)
    {
    }

    double EnergyDecay(double energy, double rnd, double density) override
    {
        auto rndd = -std::log(rnd) * density;
        auto rnddMin = decay_integral.Calculate(energy, lower_lim) / lifetime;
        if (rndd >= rnddMin)
            return lower_lim;
        return decay_integral.GetUpperLimit(energy, rndd * lifetime);
    }
};

template <typename T>
std::unique_ptr<Decay> make_decay(
    T&& cross, const ParticleDef& p_def, bool interpolate = true)
{
    if (interpolate)
        return PROPOSAL::make_unique<DecayBuilder<UtilityInterpolant, T>>(
            std::forward<T>(cross), p_def);
    return PROPOSAL::make_unique<DecayBuilder<UtilityIntegral, T>>(
            std::forward<T>(cross), p_def);
}
} // namespace PROPOSAL
