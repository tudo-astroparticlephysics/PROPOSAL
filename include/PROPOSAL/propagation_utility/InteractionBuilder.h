#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include <type_traits>

namespace PROPOSAL {
template <class I> class InteractionBuilder : public Interaction {
    I interaction_integral;

public:
    I BuildInteractionIntegral()
    {
        auto disp = std::shared_ptr<Displacement>(
                make_displacement(cross_list, false));
        auto interaction_func = [this, disp](double energy) {
            return FunctionToIntegral(*disp, energy);
        };
        I integral(interaction_func, lower_lim);
        if (typeid(I) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross_list);
            integral.BuildTables("interaction", hash, interpol_def);
        };
        return integral;
    }

    template <typename Cross>
    InteractionBuilder(Cross&& cross)
        : Interaction(std::forward<Cross>(cross))
        , interaction_integral(BuildInteractionIntegral())
    {
    }

    double EnergyInteraction(double energy, double rnd) final
    {
        assert(energy >= lower_lim);
        auto rndi = -std::log(rnd);
        auto rndiMin = interaction_integral.Calculate(energy, lower_lim);
        if (rndi >= rndiMin)
            return lower_lim;
        return interaction_integral.GetUpperLimit(energy, rndi);
    }

};

template <typename T>
std::unique_ptr<Interaction> make_interaction(T&& cross, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<InteractionBuilder<UtilityInterpolant>>(
            std::forward<T>(cross));
    return PROPOSAL::make_unique<InteractionBuilder<UtilityIntegral>>(
        std::forward<T>(cross));
}

} // namespace PROPOSAL
