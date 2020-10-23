#include "PROPOSAL/propagation_utility/Interaction.h"

using namespace PROPOSAL;

Interpolant1DBuilder::Definition Interaction::interpol_def = {1000};

double Interaction::FunctionToIntegral(Displacement& disp, double energy)
{
    auto total_rate = 0.;
    for (auto& c : cross_list) {
        total_rate += c->CalculatedNdx(energy);
    }
    if(total_rate > 0)
        return disp.FunctionToIntegral(energy) * total_rate;
    return 0;
}

double Interaction::MeanFreePath(double energy)
{
    auto total_rate = 0.;
    for (auto& c : cross_list) {
        total_rate += c->CalculatedNdx(energy);
    }
    return 1 / total_rate;
}

Interaction::loss_t Interaction::SampleLoss(double energy, std::vector<rate_t> const& rates, double rnd)
{
    auto sampled_rate = rnd * std::accumulate(rates.begin(), rates.end(), 0., [](double a, rate_t r) { return a + std::get<RATE>(r); });
    for (auto& r : rates) {
        sampled_rate -= std::get<RATE>(r);
        if (sampled_rate < 0.) {
            auto loss = std::get<CROSS>(r)->CalculateStochasticLoss(
                    std::get<COMP>(r), energy, -sampled_rate);
            return std::make_tuple( std::get<CROSS>(r)->GetInteractionType(),
                    std::get<COMP>(r), loss);
        }
    }
    throw std::logic_error("Given rate is larger than overall crosssection rate.");
}

std::vector<Interaction::rate_t> Interaction::Rates(double energy){
    auto rates = std::vector<rate_t>();
    for (auto& c : cross_list) {
        auto comp_list = c->GetTargets();
        for (auto comp : comp_list)
        {
            auto rates_comp = c->CalculatedNdx(energy, comp);
            rates.emplace_back(c, comp, rates_comp);
        }
    }
    return rates;
}
