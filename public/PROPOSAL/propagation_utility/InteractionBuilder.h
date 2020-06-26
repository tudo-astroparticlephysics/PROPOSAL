#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Interaction.h"

namespace PROPOSAL {
template <class T, class Cross> class InteractionBuilder : public Interaction {
    T interaction_integral;
    Cross crosssection_list;

    using cross_type =
        typename std::decay<Cross>::type::value_type::element_type;

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

    enum { CROSS, COMP, RATE };
    tuple<InteractionType, const Component*, double> TypeInteraction(
        double energy, double rnd) override
    {
        using rates_t = tuple<cross_type*, const Component*, double>;
        auto rates = std::vector<rates_t>();
        for (auto& c : crosssection_list) {
            auto rates_comp = c->CalculatedNdx(energy);
            for (auto& r : rates_comp)
                rates.emplace_back(c.get(), r.first, r.second);
        }
        auto sum_of_rates = std::accumulate(rates.begin(), rates.end(), 0.f,
            [](double a, rates_t r) { return a + std::get<RATE>(r); });
        sum_of_rates *= rnd;
        for (auto& r : rates) {
            sum_of_rates -= std::get<RATE>(r);
            if (sum_of_rates < 0.f) {
                auto loss = std::get<CROSS>(r)->CalculateStochasticLoss(
                    *std::get<COMP>(r), energy, -sum_of_rates);
                return std::make_tuple(std::get<CROSS>(r)->GetInteractionType(),
                    std::get<COMP>(r), loss);
            }
        }
        throw std::logic_error(
            "Given rate is larger than overall crosssection rate.");
    }
};

template <typename T>
std::unique_ptr<Interaction> make_interaction(T&& cross, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<InteractionBuilder<UtilityInterpolant, T>>(
            std::forward<T>(cross));
    return PROPOSAL::make_unique<InteractionBuilder<UtilityIntegral, T>>(
        std::forward<T>(cross));
}

} // namespace PROPOSAL
