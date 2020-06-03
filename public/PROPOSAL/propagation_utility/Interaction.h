#pragma once
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class Interaction {
protected:
    double lower_lim;

    template <typename Cross, typename Disp>
    double FunctionToIntegral(Cross&&, Disp&&, double);

public:
    Interaction(double lower_lim)
        : lower_lim(lower_lim){};
    virtual ~Interaction() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    enum { INTERACTION_TYPE, LOSS };
    virtual double EnergyInteraction(double, double) = 0;
    virtual tuple<InteractionType, double> TypeInteraction(double, double)
        = 0; // ARGS: energy, rate
};

template <class T, class Cross> class InteractionBuilder : public Interaction {
    T interaction_integral;
    Cross crosssection_list;

public:
    InteractionBuilder(Cross&&);

    T BuildInteractionIntegral(Cross&&);

    double EnergyInteraction(double, double) override;
    tuple<InteractionType, double> TypeInteraction(double, double) override;
};

template <typename Cross, typename Disp>
double Interaction::FunctionToIntegral(
    Cross&& cross, Disp&& disp, double energy)
{
    auto total_rate = 0.f;
    for (auto& c : cross)
        total_rate += c->CalculatedNdx(energy);

    return disp.FunctionToIntegral(energy) * total_rate;
}

template <class T, class Cross>
InteractionBuilder<T, Cross>::InteractionBuilder(Cross&& cross)
    : Interaction(GetLowerEnergyLim(cross))
    , interaction_integral(BuildInteractionIntegral())
    , crosssection_list(cross)
{
}

template <class T, class Cross>
T InteractionBuilder<T, Cross>::BuildInteractionIntegral(Cross&& cross)
{
    auto disp = DisplacementBuilder<UtilityIntegral, Cross>(cross);
    auto interaction_func = [this, &cross, &disp](double energy) {
        return FunctionToIntegral(cross, disp, energy);
    };
    T integral(interaction_func, GetLowerLim(cross));
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = disp.GetHash(cross);
        integral.BuildTables("interaction", hash, interpol_def);
    };
    return integral;
}

template <class T, class Cross>
double InteractionBuilder<T, Cross>::EnergyInteraction(
    double energy, double rnd)
{
    assert(energy >= lower_lim);
    auto rndi = -std::log(rnd);
    auto rndiMin = interaction_integral.Calculate(energy, lower_lim, rndi);
    if (rndi >= rndiMin)
        return lower_lim;
    return interaction_integral.GetUpperLimit(energy, rndi);
}

template <class T, class Cross>
tuple<InteractionType, double> InteractionBuilder<T, Cross>::TypeInteraction(
    double energy, double rate)
{
    for (auto& c : crosssection_list) {
        auto rates = c->CalculatedNdx(energy);
        for (auto& r : rates) {
            rate -= r;
            if (rate < 0) {
                auto loss = c->CalculateStochasticLoss(r.first, energy, -rate);
                return make_tuple(c->GetInteractionType(), loss);
            }
        }
    }
    throw std::logic_error("Given rate is larger than overall crosssection rate.");
}
}
