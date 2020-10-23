#pragma once

#include "PROPOSAL/propagation_utility/Displacement.h"
#include <tuple>

namespace PROPOSAL {
class Interaction {
protected:
    using rate_t = std::tuple<std::shared_ptr<CrossSectionBase>, std::shared_ptr<const Component>, double>;
    using loss_t = std::tuple<InteractionType, std::shared_ptr<const Component>, double>;
    using crossbase_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;

    crossbase_list_t cross_list;
    double lower_lim;

    double FunctionToIntegral(Displacement&, double);

public:
    template <typename Cross>
    Interaction(Cross const& cross)
        : cross_list(std::begin(cross), std::end(cross))
        , lower_lim(CrossSectionVector::GetLowerLim(cross))
    {}
    virtual ~Interaction() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double EnergyInteraction(double, double) = 0;

    enum { INTERACTION_TYPE, COMPONENT, FRACTIONAL_LOSS };
    std::vector<rate_t> Rates(double energy);

    enum { CROSS, COMP, RATE };
    loss_t SampleLoss(double energy, std::vector<rate_t> const& rates, double rnd);

    double MeanFreePath(double energy);
};
} // namespace PROPOSAL
