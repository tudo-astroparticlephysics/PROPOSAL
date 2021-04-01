#pragma once

#include <memory>
#include <vector>

namespace PROPOSAL {
struct CrossSectionBase;
class Component;
class Displacement;
enum class InteractionType;
}

namespace PROPOSAL {
class Interaction {

protected:
    using cross_ptr = std::shared_ptr<CrossSectionBase>;
    using crosssection_list_t = std::vector<cross_ptr>;

    std::shared_ptr<Displacement> disp;
    crosssection_list_t cross_list;
    size_t hash;


public:
    Interaction(std::shared_ptr<Displacement>, crosssection_list_t const&);
    virtual ~Interaction() = default;

    virtual double EnergyInteraction(double, double) = 0;
    virtual double EnergyIntegral(double, double) = 0;
    double FunctionToIntegral(double) const;

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
    Loss SampleLoss(double energy, std::vector<Rate> const& rates, double rnd);

    double MeanFreePath(double energy);

    auto GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
