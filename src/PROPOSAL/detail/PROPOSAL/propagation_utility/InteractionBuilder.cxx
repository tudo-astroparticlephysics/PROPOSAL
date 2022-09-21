#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

InteractionBuilder::InteractionBuilder(std::shared_ptr<Displacement> _disp,
    std::vector<cross_ptr> const& _cross, std::false_type,
    bool interpolate_meanfreepath)
    : Interaction(_disp, _cross)
    , interaction_integral(std::make_unique<UtilityIntegral>(
          [this](double E) { return FunctionToIntegral(E); },
          disp->GetLowerLim(), this->GetHash()))
{
    if (interpolate_meanfreepath)
        rate_interpolant_ = InitializeRateInterpolant();
    else
        rate_interpolant_ = nullptr;
}

InteractionBuilder::InteractionBuilder(std::shared_ptr<Displacement> _disp,
    std::vector<cross_ptr> const& _cross, std::true_type,
    bool interpolate_meanfreepath)
    : Interaction(_disp, _cross)
    , interaction_integral(std::make_unique<UtilityInterpolant>(
          [this](double E) { return FunctionToIntegral(E); },
          _disp->GetLowerLim(), this->GetHash()))
{
    interaction_integral->BuildTables("inter_",
                                      InterpolationSettings::NODES_UTILITY, false);

    if (interpolate_meanfreepath)
        rate_interpolant_ = InitializeRateInterpolant();
    else
        rate_interpolant_ = nullptr;
}

InteractionBuilder::interpolant_ptr InteractionBuilder::InitializeRateInterpolant() {
    auto energy_lim = AxisBuilderDNDX::energy_limits();
    energy_lim.low = disp->GetLowerLim();
    energy_lim.up = InterpolationSettings::UPPER_ENERGY_LIM;
    energy_lim.nodes = InterpolationSettings::NODES_RATE_INTERPOLANT;
    auto energy_lim_refined = AxisBuilderDNDX::refine_definition_range(
            energy_lim, [&](double E) { return calculate_total_rate(E); });
    auto def = cubic_splines::CubicSplines<double>::Definition();
    def.f = [&](double energy) {
        return calculate_total_rate(energy);
    };
    auto axis = AxisBuilderDNDX::Create(energy_lim_refined);
    def.axis = std::move(axis);
    rate_lower_energy_lim = def.axis->GetLow();

    auto rate_interpolant_hash = this->GetHash();
    hash_combine(rate_interpolant_hash,
                 InterpolationSettings::NODES_RATE_INTERPOLANT,
                 InterpolationSettings::UPPER_ENERGY_LIM);

    return std::make_shared<interpolant_t>(
            std::move(def), std::string(InterpolationSettings::TABLES_PATH),
            std::string("rates_") + std::to_string(rate_interpolant_hash) + std::string(".dat"));
}

double InteractionBuilder::EnergyInteraction(double energy, double rnd)
{
    assert(energy >= disp->GetLowerLim());
    auto rndi = -std::log(rnd);
    auto rndiMin = interaction_integral->Calculate(energy, disp->GetLowerLim());
    if (rndi >= rndiMin)
        return disp->GetLowerLim();
    return interaction_integral->GetUpperLimit(energy, rndi);
}

double InteractionBuilder::EnergyIntegral(double E_i, double E_f) {
    return interaction_integral->Calculate(E_i, E_f);
}

double InteractionBuilder::MeanFreePath(double energy) {
    if (rate_interpolant_) {
        if (energy < rate_lower_energy_lim)
            return INF;
        auto rate = rate_interpolant_->evaluate(energy);
        if (rate < 0) {
            Logging::Get("proposal.interaction")->warn(
                    "Negative MeanFreePath detected at energy {} MeV. Returning INF instead.", energy);
            return INF;
        }
        return 1. / rate;
    }
    return 1. / calculate_total_rate(energy);
}

namespace PROPOSAL {
std::unique_ptr<Interaction> make_interaction(
    std::shared_ptr<Displacement> disp,
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross,
    bool interpolate_interaction_integral, bool interpolate_meanfreepath)
{
    auto inter = std::unique_ptr<Interaction>();
    if (interpolate_interaction_integral)
        inter = std::make_unique<InteractionBuilder>(
                disp, cross, std::true_type {}, interpolate_meanfreepath);
    else
        inter = std::make_unique<InteractionBuilder>(
                disp, cross, std::false_type {}, interpolate_meanfreepath);
    return inter;
}

std::unique_ptr<Interaction> make_interaction(
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross,
    bool interpolate_interaction_integral, bool interpolate_meanfreepath)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    return make_interaction(disp, cross, interpolate_interaction_integral,
                            interpolate_meanfreepath);
}
} // namespace PROPOSAL
