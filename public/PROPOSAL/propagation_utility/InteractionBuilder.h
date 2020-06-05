#include "PROPOSAL/propagation_utility/Interaction.h"

namespace PROPOSAL {
template <class T, class Cross> class InteractionBuilder : public Interaction {
    T interaction_integral;
    Cross crosssection_list;

public:
    InteractionBuilder(Cross&& cross)
        : Interaction(CrossSectionVector::GetLowerLim(cross))
        , interaction_integral(BuildInteractionIntegral(cross))
        , crosssection_list(cross)
    {
    }

    T BuildInteractionIntegral(Cross&& cross)
    {
        auto disp = DisplacementBuilder<UtilityIntegral, Cross>(cross);
        auto interaction_func = [this, &cross, &disp](double energy) {
            return FunctionToIntegral(cross, disp, energy);
        };
        T integral(interaction_func, CrossSectionVector::GetLowerLim(cross));
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross);
            integral.BuildTables("interaction", hash, interpol_def);
        };
        return integral;
    }

    double EnergyInteraction(double energy, double rnd) override
    {
        assert(energy >= lower_lim);
        auto rndi = -std::log(rnd);
        auto rndiMin = interaction_integral.Calculate(energy, lower_lim);
        if (rndi >= rndiMin)
            return lower_lim;
        return interaction_integral.GetUpperLimit(energy, rndi);
    }

    tuple<InteractionType, double> TypeInteraction(
        double energy, double rate) override
    {
        for (auto& c : crosssection_list) {
            auto rates = c->CalculatedNdx(energy);
            for (auto& r : rates) {
                rate -= r.second;
                if (rate < 0) {
                    auto loss
                        = c->CalculateStochasticLoss(*r.first, energy, -rate);
                    return std::make_tuple(c->GetInteractionType(), loss);
                }
            }
        }
        throw std::logic_error(
            "Given rate is larger than overall crosssection rate.");
    }
};

template <typename T>
std::unique_ptr<Interaction> make_Interaction(T&& cross, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<
            InteractionBuilder<UtilityInterpolant, T>>(std::forward<T>(cross));
    return PROPOSAL::make_unique<InteractionBuilder<UtilityIntegral, T>>(
        std::forward<T>(cross));
}

} // namespace PROPOSAL
