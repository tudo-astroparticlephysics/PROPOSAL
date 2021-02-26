#pragma once
#include "PROPOSAL/propagation_utility/Interaction.h"

#include <type_traits>
#include <vector>

namespace PROPOSAL {
class UtilityIntegral;
} // namespace PROPOSAL

namespace PROPOSAL {
class InteractionBuilder : public Interaction {
    std::unique_ptr<UtilityIntegral> interaction_integral;

public:
    InteractionBuilder(std::shared_ptr<Displacement>,
        crosssection_list_t const&, std::false_type);

    InteractionBuilder(std::shared_ptr<Displacement>,
        crosssection_list_t const&, std::true_type);

    double EnergyInteraction(double energy, double rnd) final;
    double EnergyIntegral(double E_i, double E_f) final;
};
} // namespace PROPOSAL

namespace PROPOSAL {
std::unique_ptr<Interaction> make_interaction(std::shared_ptr<Displacement>,
    std::vector<std::shared_ptr<CrossSectionBase>> const&, bool);

std::unique_ptr<Interaction> make_interaction(
    std::vector<std::shared_ptr<CrossSectionBase>> const&, bool);
} // namespace PROPOSAL
