#pragma once
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/Interaction.h"

namespace PROPOSAL {
template <class I> class InteractionBuilder : public Interaction {
    I interaction_integral;

    void build_tables() {};

public:

    template <typename Cross>
    InteractionBuilder(std::shared_ptr<Displacement> _disp, Cross&& _cross)
        : Interaction(_disp, std::forward<Cross>(_cross))
        , interaction_integral(
              [this](double E) { return FunctionToIntegral(E); },
              disp->GetLowerLim(), this->GetHash())
    {
    }

    double EnergyInteraction(double energy, double rnd) final
    {
        assert(energy >= disp->GetLowerLim());
        auto rndi = -std::log(rnd);
        auto rndiMin
            = interaction_integral.Calculate(energy, disp->GetLowerLim());
        if (rndi >= rndiMin)
            return disp->GetLowerLim();
        return interaction_integral.GetUpperLimit(energy, rndi);
    }
};

template <typename T>
auto make_interaction(
    std::shared_ptr<Displacement> disp, T&& cross, bool interpolate)
{
    auto inter = std::unique_ptr<Interaction>();
    if (interpolate)
        inter.reset(new InteractionBuilder<UtilityInterpolant>(
            disp, std::forward<T>(cross)));
    else
        inter.reset(new InteractionBuilder<UtilityIntegral>(
            disp, std::forward<T>(cross)));
    return inter;
}

template <typename T> auto make_interaction(T&& cross, bool interpolate)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    return make_interaction(disp, cross, interpolate);
}

} // namespace PROPOSAL
