#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <iostream>

using namespace PROPOSAL;

InteractionBuilder::InteractionBuilder(std::shared_ptr<Displacement> _disp,
    std::vector<cross_ptr> const& _cross, std::false_type)
    : Interaction(_disp, _cross)
    , interaction_integral(std::make_unique<UtilityIntegral>(
          [this](double E) { return FunctionToIntegral(E); }, GetLowerLim(),
          this->GetHash()))
{
}

InteractionBuilder::InteractionBuilder(std::shared_ptr<Displacement> _disp,
    std::vector<cross_ptr> const& _cross, std::true_type)
    : Interaction(_disp, _cross)
    , interaction_integral(std::make_unique<UtilityInterpolant>(
          [this](double E) { return FunctionToIntegral(E); }, GetLowerLim(),
          this->GetHash()))
{
    interaction_integral->BuildTables("inter_", 500, false);
}

double InteractionBuilder::EnergyInteraction(double energy, double rnd)
{
    /* assert(energy >= disp->GetLowerLim()); */
    auto rndi = -std::log(rnd);
    auto rndiMin = EnergyIntegral(energy, GetLowerLim());

    if (rndi >= rndiMin)
        return GetLowerLim();

    return interaction_integral->GetUpperLimit(energy, rndi);
}

double InteractionBuilder::EnergyIntegral(double E_i, double E_f)
{
    if (E_f < GetLowerLim())
        return interaction_integral->Calculate(E_i, GetLowerLim());

    return interaction_integral->Calculate(E_i, E_f);
}

namespace PROPOSAL {
std::unique_ptr<Interaction> make_interaction(
    std::shared_ptr<Displacement> disp,
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross, bool interpol)
{
    auto inter = std::unique_ptr<Interaction>();
    if (interpol)
        inter = std::make_unique<InteractionBuilder>(
            disp, cross, std::true_type {});
    else
        inter = std::make_unique<InteractionBuilder>(
            disp, cross, std::false_type {});
    return inter;
}

std::unique_ptr<Interaction> make_interaction(
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross, bool interpol)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    return make_interaction(disp, cross, interpol);
}
} // namespace PROPOSAL
