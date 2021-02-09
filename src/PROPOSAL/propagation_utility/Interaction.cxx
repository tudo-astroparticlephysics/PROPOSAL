#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

Interaction::Interaction(
    std::shared_ptr<Displacement> _disp, std::vector<cross_ptr> const& _cross)
    : disp(_disp)
    , cross_list(_cross)
{
}

double Interaction::FunctionToIntegral(double energy) const
{
    auto total_rate = 0.;
    for (auto& c : cross_list)
        total_rate += c->CalculatedNdx(energy);
    if (total_rate > 0)
        return disp->FunctionToIntegral(energy) * total_rate;
    return 0;
}

double Interaction::MeanFreePath(double energy)
{
    auto total_rate = 0.;
    for (auto& c : cross_list)
        total_rate += c->CalculatedNdx(energy);
    return 1 / total_rate;
}

Interaction::Loss Interaction::SampleLoss(
    double energy, std::vector<Rate> const& rates, double rnd)
{
    auto sampled_rate = rnd
        * std::accumulate(rates.begin(), rates.end(), 0.,
            [](double a, Rate r) { return a + r.rate; });
    for (auto& r : rates) {
        sampled_rate -= r.rate;
        if (sampled_rate < 0.) {
            auto loss = r.crosssection->CalculateStochasticLoss(
                r.comp_hash, energy, -sampled_rate);
            return { r.crosssection->GetInteractionType(), r.comp_hash, loss };
        }
    }
    throw std::logic_error(
        "Given rate is larger than overall crosssection rate.");
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
