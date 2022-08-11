#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

#include <sstream>
#include <numeric>

using namespace PROPOSAL;

Interaction::Interaction(
    std::shared_ptr<Displacement> _disp, std::vector<cross_ptr> const& _cross)
    : disp(_disp)
    , cross_list(_cross)
    , hash(CrossSectionVector::GetHash(cross_list))
{
    if (cross_list.size() < 1)
        throw std::invalid_argument("At least one crosssection is required.");
}

double Interaction::FunctionToIntegral(double energy) const
{
    auto total_rate = calculate_total_rate(energy);
    if (total_rate > 0)
        return disp->FunctionToIntegral(energy) * total_rate;
    return 0;
}

double Interaction::calculate_total_rate(double energy) const {
    auto total_rate = 0.;
    for (auto& c : cross_list)
        total_rate += c->CalculatedNdx(energy);
    return total_rate;
}

Interaction::Loss Interaction::SampleLoss(
    double energy, std::vector<Rate> const& rates, double rnd)
{
    auto overall_rate = std::accumulate(rates.begin(), rates.end(), 0.,
        [](double a, Rate r) { return a + r.rate; });
    auto sampled_rate = rnd * overall_rate;
    for (auto& r : rates) {
        sampled_rate -= r.rate;
        if (sampled_rate < 0.) {
            auto loss = r.crosssection->CalculateStochasticLoss(
                r.comp_hash, energy, -sampled_rate);
            return { r.crosssection->GetInteractionType(), r.comp_hash, loss };
        }
    }

    if (overall_rate == 0.) {
        Logging::Get("proposal.interaction")->warn(
                "No stochastic interaction possible for initial energy {} MeV.",
                energy);
        return {InteractionType::Undefined, 0, 0};
    }

    std::stringstream ss;
    ss << "Given rate (" << std::to_string(sampled_rate)
       << ") for given energy (" << std::to_string(energy)
       << ") by drawn random number (" << std::to_string(rnd)
       << ") is larger than the overall crosssection rate ("
       << std::to_string(overall_rate) << ").";

    throw std::logic_error(ss.str());
}

std::vector<Interaction::Rate> Interaction::Rates(double energy)
{
    auto rates = std::vector<Rate>();
    for (auto& c : cross_list) {
        auto dndx_per_target = c->CalculatedNdx_PerTarget(energy);
        for (auto dndx : dndx_per_target)
            rates.emplace_back(
                Interaction::Rate { c, dndx.first, dndx.second });
    }
    return rates;
}
