/* #include "PROPOSAL/propagation_utility/Interaction.h" */

/* using namespace PROPOSAL; */
/* using std::make_tuple; */

/* namespace PROPOSAL { */
/* Interpolant1DBuilder::Definition interaction_interpol_def; */
/* } // namespace PROPOSAL */

/* Interaction::Interaction(CrossSectionList cross) */
/*     : cross(cross) */
/*     , displacement(cross) */
/*     , lower_lim(std::numeric_limits<double>::max()) */
/* { */
/*     if (cross.size() < 1) */
/*         throw std::invalid_argument("at least one crosssection is required."); */

/*     for (auto c : cross) */
/*         lower_lim = std::min(lower_lim, c->GetLowerEnergyLimit()); */
/* } */

/* unordered_map<InteractionType, comp_rates> Interaction::CalculateRates(double energy) */
/* { */
/*     unordered_map<InteractionType, comp_rates> rates; */
/*     for (const auto& c: cross) */
/*         rates[c->parametrization_->GetType()] = c->CalculatedNdx(energy); */

/*     return rates; */
/* } */

/* double Interaction::FunctionToIntegral(double energy) */
/* { */
/*     auto interaction_rates = CalculateRates(energy); */

/*     auto total_rate = (double)0; */
/*     for (const auto& comp_rates : interaction_rates) { */
/*         for (const auto& rate : comp_rates.second) */
/*             total_rate += rate.second; */
/*     } */

/*     return displacement.FunctionToIntegral(energy) * total_rate; */
/* } */

/* tuple<double, InteractionType> Interaction::SampleStochasticLoss( */
/*     double energy, double rnd) */
/* { */
/*     auto interaction_rates = CalculateRates(energy); */
/*     auto total_rate = (double)0; */
/*     for (const auto& comp_rates : interaction_rates) { */
/*         for (const auto& rate : comp_rates.second) */
/*             total_rate += rate.second; */
/*     } */

/*     auto sampled_rate = total_rate * rnd; */
/*     auto inter = begin(cross); */
/*     for (const auto& inter_rate : interaction_rates) { */
/*         for (const auto& rate : inter_rate.second) { */
/*             sampled_rate -= rate.second; */
/*             if (sampled_rate < 0) { */
/*                 auto loss = (*inter)->CalculateStochasticLoss( */
/*                     energy, -sampled_rate, rate.first); */
/*                 return make_tuple(loss, inter_rate.first); */
/*             } */
/*         } */
/*     } */

/*     throw std::logic_error("Sampling interaction rate goes wrong. Maybe the " */
/*                            "random number is not in range [0,1] "); */
/* } */
