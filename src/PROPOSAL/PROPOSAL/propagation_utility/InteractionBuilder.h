#pragma once
#include "PROPOSAL/propagation_utility/Interaction.h"
#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

#include <type_traits>
#include <vector>

namespace PROPOSAL {
class UtilityIntegral;
} // namespace PROPOSAL

namespace PROPOSAL {
class InteractionBuilder : public Interaction {
    std::unique_ptr<UtilityIntegral> interaction_integral;

    using interpolant_t
        = cubic_splines::Interpolant<cubic_splines::CubicSplines<double>>;
    using interpolant_ptr = std::shared_ptr<interpolant_t>;
    interpolant_ptr rate_interpolant_;
    double rate_lower_energy_lim;

    interpolant_ptr InitializeRateInterpolant();

public:
    InteractionBuilder(std::shared_ptr<Displacement>,
        crosssection_list_t const&, std::false_type, bool);

    InteractionBuilder(std::shared_ptr<Displacement>,
        crosssection_list_t const&, std::true_type, bool);

    double EnergyInteraction(double energy, double rnd) final;
    double EnergyIntegral(double E_i, double E_f) final;

    double MeanFreePath(double energy) final;

};
} // namespace PROPOSAL

namespace PROPOSAL {
std::unique_ptr<Interaction> make_interaction(std::shared_ptr<Displacement>,
    std::vector<std::shared_ptr<CrossSectionBase>> const&, bool, bool = false);

std::unique_ptr<Interaction> make_interaction(
    std::vector<std::shared_ptr<CrossSectionBase>> const&, bool, bool = false);
} // namespace PROPOSAL
