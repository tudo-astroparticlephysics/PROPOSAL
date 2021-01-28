#pragma once

#include <tuple>
#include <memory>
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

namespace PROPOSAL {
    class CrossSectionBase;
    namespace Components {
        class Component;
    }
    class Displacement;
}

namespace PROPOSAL {
class Interaction {
    size_t hash;

protected:
    using comp_ptr = std::shared_ptr<const Components::Component>;
    using cross_ptr = std::shared_ptr<CrossSectionBase>;
    using rate_t = std::tuple<cross_ptr, comp_ptr, double>;
    using loss_t = std::tuple<InteractionType, comp_ptr, double>;

    std::shared_ptr<Displacement> disp;
    std::vector<cross_ptr> cross_list;

    double FunctionToIntegral(double);

public:
    template <typename Cross>
    Interaction(std::shared_ptr<Displacement> _disp, Cross const& _cross)
        : disp(_disp)
        , cross_list(std::begin(_cross), std::end(_cross))
    {
    }
    virtual ~Interaction() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double EnergyInteraction(double, double) = 0;

    enum { INTERACTION_TYPE, COMPONENT, FRACTIONAL_LOSS };
    std::vector<rate_t> Rates(double energy);

    enum { CROSS, COMP, RATE };
    loss_t SampleLoss(
        double energy, std::vector<rate_t> const& rates, double rnd);

    double MeanFreePath(double energy);
    auto GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
