
#include "PROPOSAL/Logging.h"
#include <PROPOSAL/crossection/factories/PhotoPairFactory.h>

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;

/******************************************************************************
 *                            Propagation utility                              *
 ******************************************************************************/

Utility::Definition::Definition(std::shared_ptr<Crosssections> corsssection,
    std::shared_ptr<Scattering> scattering = nullptr,
    std::shared_ptr<InterpolationDef> inter_def = nullptr)
    : crosssection(crosssection)
    , scattering(scattering)
    , inter_def(inter_def)
{
    if (!inter_def) {
        log_warn("No interpolation definition defined. Integral will not be "
                 "approximate by interpolants. Performance will be poor.");

        displacement_calculator.reset(
            new UtilityIntegralDisplacement(crosssection));
        interaction_calculator.reset(
            new UtilityIntegralInteraction(crosssection));
        decay_calculator.reset(new UtilityIntegralDecay(crosssection))
    } else {
        displacement_calculator.reset(
            new UtilityInterpolantDisplacement(crosssection, inter_def));
        interaction_calculator.reset(
            new UtilityInterpolantInteraction(crosssection, inter_def));
        decay_calculator.reset(new UtilityInterpolantDecay(crosssection, inter_def))
    }

    if (!scattering) {
        log_debug("No scattering defined. Partilce will only be deflected if a "
                 "stochastic deflection is in the crosssection implemented.");
    }

    if (!cont_rand) {
        log_debug("No continous randomization used.");
    } else {
        if (!inter_def) {
            cont_rand.reset(new UtilityIntegralContRand(crosssection))
        } else {
            cont_rand.reset(new UtilityInterpolantContRand(crosssection, inter_def))
        }
    }

    if (!exact_time) {
        log_debug("No exact time calculator used.");
    }
    } else {
        if (!inter_def) {
            cont_rand.reset(new UtilityIntegralTime(crosssection))
        } else {
            cont_rand.reset(new UtilityInterpolantTime(crosssection, inter_def))
        }
    }
}

std::ostream& PROPOSAL::operator<<(
    std::ostream& os, PROPOSAL::Utility::Definition const& util_definition)
{
    std::stringstream ss;
    ss << " Utility Definition (" << &util_definition << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    for (const auto& crosssection : util_definition->crosssection) {
        os << crosssection << std::endl;
    }
    if (scattering) {
        os << scattering << std::endl;
    };
    if (inter_def) {
        os << inter_def << std::endl;
    }
    if (cont_rand) {
        os << cont_rand << std::endl;
    }
    if (exact_time) {
        os << exact_time << std::endl;
    } ;
    os << Helper::Centered(60, "");
    return os;
}

// -------------------------------------------------------------------------
// // Constructors
// -------------------------------------------------------------------------

Utility::Utility(std::unique_ptr<Utility::Definition> utility_def)
    : utility_def(std::move(utility_def))
{
}

std::pair<double, int> Utility::StochasticLoss(
    double particle_energy, double rnd1, double rnd2, double rnd3)
{
    double total_rate = 0;
    double total_rate_weighted = 0;
    double rates_sum = 0;
    std::vector<double> rates;
    rates.resize(crosssections_.size());

    // return 0 and unknown, if there is no interaction
    std::pair<double, int> energy_loss;
    energy_loss.first = 0.;
    energy_loss.second = 0;

    for (unsigned int i = 0; i < crosssections_.size(); i++) {
        rates[i] = crosssections_[i]->CalculatedNdx(particle_energy, rnd2);
        total_rate += rates[i];
    }

    total_rate_weighted = total_rate * rnd1;

    log_debug("Total rate = %f, total rate weighted = %f", total_rate,
        total_rate_weighted);

    for (unsigned int i = 0; i < rates.size(); i++) {
        rates_sum += rates[i];

        if (rates_sum >= total_rate_weighted) {
            energy_loss.first = crosssections_[i]->CalculateStochasticLoss(
                particle_energy, rnd2, rnd3);
            energy_loss.second = crosssections_[i]->GetTypeId();

            break;
        }
    }

    return energy_loss;
}
