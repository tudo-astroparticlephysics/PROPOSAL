#include "PROPOSAL/propagation_utility/Interaction.h"

using namespace PROPOSAL;

namespace PROPOSAL {
Interpolant1DBuilder::Definition interaction_interpol_def;
} // namespace PROPOSAL

Interaction::Interaction(CrossSectionList cross)
    : cross(cross)
    , lower_lim(std::numeric_limits<double>::max())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");

    for (auto c : cross)
        lower_lim = std::min(lower_lim, c->GetLowerEnergyLimit());
}

std::shared_ptr<CrossSection> Interaction::TypeInteraction(
    double energy, const std::array<double, 2>& rnd)
{
    /* assert(energy >= lower_lim); */
    /* std::vector<double> rates; */
    /* for (const auto& c : cross) */
    /*     rates.push_back(c->CalculatedNdx(energy, rnd[1])); */

    /* auto total_rate = std::accumulate(rates.begin(), rates.end(), 0.0); */
    /* log_debug("Total rate = %f, total rate weighted = %f", total_rate, */
    /*     total_rate * rnd[0]); */

    /* double rates_sum = 0.; */
    /* for (size_t i = 0; i < rates.size(); i++) { */
    /*     rates_sum += rates[i]; */
    /*     if (rates_sum >= total_rate * rnd[0]) */
    /*         return cross.at(i); */
    /* } */

    throw std::logic_error(
        "Something went wrong during the total rate calculation.");
}
