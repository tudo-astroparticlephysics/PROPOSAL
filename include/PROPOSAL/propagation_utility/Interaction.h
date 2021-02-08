#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/particle/Particle.h"
#include <memory>
#include <tuple>

namespace PROPOSAL {
struct CrossSectionBase;
class Component;
class Displacement;
}

namespace PROPOSAL {
class Interaction {
    size_t hash;

protected:
    using comp_ptr = std::shared_ptr<const Component>;
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

    struct Rate {
        cross_ptr crosssection;
        size_t comp_hash;
        double rate;
    };
    std::vector<Rate> Rates(double energy);

    struct Loss {
        InteractionType type;
        size_t comp_hash;
        double v_loss;
    };
    Loss SampleLoss(
        double energy, std::vector<Rate> const& rates, double rnd);

    double MeanFreePath(double energy);
    auto GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
