#include "PROPOSAL/propagation_utility/Interaction.h"

using namespace PROPOSAL;

namespace PROPOSAL {
Interpolant1DBuilder::Definition interaction_interpol_def;
} // namespace PROPOSAL

Interaction::Interaction(CrossSectionList cross)
    : cross(cross)
    , displacement(cross)
    , lower_lim(std::numeric_limits<double>::max())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");

    for (auto c : cross)
        lower_lim = std::min(lower_lim, c->GetLowerEnergyLimit());
}

vector<vector<double>> Interaction::CalculateRates(double energy)
{
    vector<vector<double>> rates;
    for (const auto& crosssection : cross)
        rates.emplace_back(crosssection->CalculatedNdx(energy));
    return rates;
}

double Interaction::FunctionToIntegral(double energy)
{
    auto interaction_rates = CalculateRates(energy);

    auto total_rate = (double)0;
    for (const auto& comp_rates : interaction_rates)
        total_rate += accumulate(comp_rates.begin(), comp_rates.end(), 0);

    return displacement.FunctionToIntegral(energy) * total_rate;
}

tuple<double, InteractionType> Interaction::SampleStochasticLoss(double energy, double rnd)
{
    auto interaction_rates = CalculateRates(energy);

    auto total_rate = (double)0;
    for (const auto& comp_rates : interaction_rates)
        total_rate += accumulate(comp_rates.begin(), comp_rates.end(), 0);

    auto sampled_rate = total_rate * rnd;

    auto interaction = cross.begin()
    for (const auto& comp_rates : interaction_rates){
        auto comp = interaction->GetComponents();
        for (const auto & rate : comp_rates) {
            sampled_rate -= rate;

            if(sampled_rate < 0){

            }
        }
        ++interaction;
    }

}
