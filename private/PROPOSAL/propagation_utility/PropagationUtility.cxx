
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"
#include <PROPOSAL/crossection/factories/PhotoPairFactory.h>

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;
using std::tuple;

/******************************************************************************
 *                            Propagation utility                              *
 ******************************************************************************/

Utility::Definition::Definition(CrossSectionList cross,
    const ParticleDef& p_def, std::shared_ptr<Scattering> scattering = nullptr,
    std::shared_ptr<InterpolationDef> interpol_def = nullptr)
    : cross(cross)
    , scattering(scattering)
    , inter_def(interpol_def)
{
    if (!interpol_def) {
        log_warn("No interpolation definition defined. Integral will not be "
                 "approximate by interpolants. Performance will be poor.");

        displacement_calc.reset(new UtilityIntegralDisplacement(cross, p_def));
        interaction_calc.reset(new UtilityIntegralInteraction(cross, p_def));
        decay_calc.reset(new UtilityIntegralDecay(cross, p_def));
    } else {
        displacement_calc.reset(
            new UtilityInterpolantDisplacement(cross, p_def, interpol_def));
        interaction_calc.reset(
            new UtilityInterpolantInteraction(cross, p_def, interpol_def));
        decay_calc.reset(
            new UtilityInterpolantDecay(cross, p_def, interpol_def));
    }

    if (!scattering) {
        log_debug("No scattering defined. Particle will only be deflected in "
                  "stochastic interactions");
    }

    if (!cont_rand) {
        log_debug("No continuous randomization used.");
    } else {
        if (!inter_def) {
            cont_rand.reset(new UtilityIntegralContRand(cross, p_def));
        } else {
            cont_rand.reset(
                new UtilityInterpolantContRand(cross, p_def, interpol_def));
        }
    }

    if (!exact_time) {
        log_debug("Exact time calculation disabled: Velocity of particles will "
                  "be approximated as speed of light");
    } else {
        if (!inter_def) {
            cont_rand.reset(new UtilityIntegralTime(cross, p_def));
        } else {
            cont_rand.reset(
                new UtilityInterpolantTime(cross, p_def, interpol_def));
        }
    }
}

/* std::ostream& PROPOSAL::operator<<( */
/*     std::ostream& os, PROPOSAL::Utility::Definition const& util_definition)
 */
/* { */
/*     std::stringstream ss; */
/*     ss << " Utility Definition (" << &util_definition << ") "; */
/*     os << Helper::Centered(60, ss.str()) << '\n'; */

/*     for (const auto& crosssection : util_definition.cross) { */
/*         os << crosssection << std::endl; */
/*     } */
/*     if (util_definition.scattering) { */
/*         os << util_definition.scattering << std::endl; */
/*     }; */
/*     if (util_definition.inter_def) { */
/*         os << util_definition.inter_def << std::endl; */
/*     } */
/*     if (util_definition.cont_rand) { */
/*         os << util_definition.cont_rand << std::endl; */
/*     } */
/*     if (util_definition.exact_time) { */
/*         os << util_definition.exact_time << std::endl; */
/*     }; */
/*     os << Helper::Centered(60, ""); */
/*     return os; */
/* } */

// -------------------------------------------------------------------------
// // Constructors
// -------------------------------------------------------------------------

Utility::Utility(std::unique_ptr<Utility::Definition> utility_def)
    : utility_def(std::move(utility_def))
{
}

std::shared_ptr<CrossSection> Utility::TypeInteraction(
    double energy, const std::array<double, 2>& rnd)
{
    std::vector<double> rates;
    for (const auto& crosssection : utility_def->cross)
        rates.push_back(crosssection->CalculatedNdx(energy, rnd[1]));

    double total_rate{ std::accumulate(rates.begin(), rates.end(), 0.0) };
    log_debug("Total rate = %f, total rate weighted = %f", total_rate,
        total_rate * rnd[0]);

    double rates_sum = 0;
    for (size_t i = 0; i < rates.size(); i++) {
        rates_sum += rates[i];
        if (rates_sum >= total_rate * rnd[0])
            return utility_def->cross.at(i);
    }

    throw std::logic_error(
        "Something went wrong during the total rate calculation.");
}

double Utility::EnergyStochasticloss(
    CrossSection& crosssection, double energy, const std::array<double, 2>& rnd)
{
    auto aux = crosssection.CalculateStochasticLoss(energy, rnd[0], rnd[1]);
    return aux;
}

double Utility::EnergyDecay(double energy, double rnd)
{
    auto rndd = -std::log(rnd);
    auto rnddMin = 0;

    rnddMin = utility_def->decay_calc->Calculate(energy, utility_def->decay_calc->GetMass(), rndd);

    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0)
        return utility_def->decay_calc->GetMass();

    return utility_def->decay_calc->GetUpperLimit(energy, rndd);
}

double Utility::EnergyInteraction(double energy, double rnd)
{
    auto rndi = -std::log(rnd);
    auto rndiMin = 0.;

    // solving the tracking integral
    rndiMin = utility_def->interaction_calc->Calculate(energy, utility_def->interaction_calc->GetMass(), rndi);

    if (rndi >= rndiMin || rndiMin <= 0)
        return utility_def->interaction_calc->GetMass();

    return utility_def->interaction_calc->GetUpperLimit(energy, rndi);
}

double Utility::EnergyRandomize(
    double initial_energy, double final_energy, double rnd)
{
    if (utility_def->cont_rand) {
        auto variance = utility_def->cont_rand->Calculate(initial_energy, final_energy, 0.0);
        return SampleFromGaussian(final_energy, variance, rnd, utility_def->cont_rand->GetMass(), initial_energy);
    } else {
        return final_energy;
    }
}

double Utility::TimeElapsed(double initial_energy, double final_energy, double distance)
{
    if (utility_def->exact_time) {
        return utility_def->exact_time->Calculate(
            initial_energy, final_energy, distance);
    } else {
        return distance / SPEED;
    }
}

tuple<Vector3D, Vector3D> Utility::DirectionsScatter(double displacement,
    double initial_energy, double final_energy, const Vector3D& position,
    const Vector3D& direction, const std::array<double, 4>& rnd)
{
    return utility_def->scattering->Scatter(
        displacement, initial_energy, final_energy, position, direction, rnd);
}

std::pair<double, double> DirectionDeflect(CrossSection& crosssection, double particle_energy, double loss_energy)
{
    return crosssection.StochasticDeflection(particle_energy, loss_energy);
}

double Utility::LengthContinuous(
    double initial_energy, double final_energy, double border_length)
{
    return utility_def->displacement_calc->Calculate(
        initial_energy, final_energy, border_length);
}
