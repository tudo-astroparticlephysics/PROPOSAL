#pragma once

#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include <tuple>

namespace PROPOSAL {
class Interaction {
protected:
    double lower_lim;

    template <typename Cross, typename Disp>
    double FunctionToIntegral(Cross&& cross, Disp&& disp, double energy)
    {
        auto total_rate = 0.f;
        for (auto& c : cross) {
            auto rates = c->CalculatedNdx(energy);
            for (auto& r : rates)
                total_rate += r.second;
        }
        if(total_rate > 0)
            return disp.FunctionToIntegral(cross, energy) * total_rate;
        return 0;
    }

    template <typename Cross>
    double CalculateMeanFreePath(Cross&& cross, double energy)
    {
        auto total_rate = 0.f;
        for (auto& c : cross) {
            auto rates = c->CalculatedNdx(energy);
            for (auto& r : rates)
                total_rate += r.second;
        }
        return 1 / total_rate;
    }

public:
    Interaction(double lower_lim)
        : lower_lim(lower_lim){};
    virtual ~Interaction() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    enum { INTERACTION_TYPE, LOSS };
    virtual double EnergyInteraction(double, double) = 0;
    virtual tuple<InteractionType, std::shared_ptr<const Component>, double> TypeInteraction(
        double, double)
        = 0; // ARGS: energy, rate
    virtual double MeanFreePath(double) = 0;
};
} // namespace PROPOSAL
