
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"

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

// Utility::Definition::Definition(CrossSectionList cross,
//     const ParticleDef& p_def, std::shared_ptr<Scattering> scattering = nullptr,
//     std::shared_ptr<InterpolationDef> inter_def = nullptr)
//     : scattering(scattering)
// {
//     if (!inter_def) {
//         log_warn("No interpolation definition defined. Integral will not be "
//                  "approximate by interpolants. Performance will be poor.");
//
//         displacement_calc.reset(new DisplacementBuilder<UtilityIntegral>(cross));
//         interaction_calc.reset(new InteractionBuilder<UtilityIntegral>(cross));
//         decay_calc.reset(new DecayBuilder<UtilityIntegral>(cross, p_def));
//     } else {
//         displacement_calc.reset(new DisplacementBuilder<UtilityInterpolant>(cross));
//         interaction_calc.reset(new InteractionBuilder<UtilityInterpolant>(cross));
//         decay_calc.reset(new DecayBuilder<UtilityInterpolant>(cross, p_def));
//     }
//
//     if (!scattering) {
//         log_debug("No scattering defined. Particle will only be deflected in "
//                   "stochastic interactions");
//     }
//
//     if (!cont_rand) {
//         log_debug("No continuous randomization used.");
//     } else {
//         if (!inter_def) {
//             cont_rand.reset(new ContRandBuilder<UtilityIntegral>(cross));
//         } else {
//             cont_rand.reset(new ContRandBuilder<UtilityInterpolant>(cross));
//         }
//     }
//
//     if (!time_calc) {
//         log_debug("Exact time calculation disabled: Velocity of particles will "
//                   "be approximated as speed of light");
//     } else {
//         if (!inter_def) {
//             time_calc.reset(new ExactTimeBuilder<UtilityIntegral>(cross, p_def));
//         } else {
//             time_calc.reset(new ExactTimeBuilder<UtilityInterpolant>(cross, p_def));
//         }
//     }
// }

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

// double Utility::EnergyStochasticloss(
//     CrossSection& crosssection, double energy, const std::array<double, 2>& rnd)
// {
//     throw std::logic_error("Not implemented jet!");
    /* auto aux = crosssection.CalculateStochasticLoss(energy, rnd[0], rnd[1]); */
    /* return aux; */
// }
//
// double Utility::EnergyDecay(double energy, double rnd)
// {
//     return utility_def->decay_calc->EnergyDecay(energy, rnd);
// }
//
// double Utility::EnergyInteraction(double energy, double rnd)
// {
//     return utility_def->interaction_calc->UpperLimitEnergyIntegral(energy, rnd);
// }
//
// double Utility::EnergyRandomize(
//     double initial_energy, double final_energy, double rnd)
// {
//     if (utility_def->cont_rand) {
//         final_energy = utility_def->cont_rand->EnergyRandomize(initial_energy, final_energy, rnd);
//     }
//     return final_energy;
// }
//
// double Utility::TimeElapsed(double initial_energy, double final_energy)
// {
//     return utility_def->time_calc->TimeElapsed(initial_energy, final_energy);
// }
//
// tuple<Vector3D, Vector3D> Utility::DirectionsScatter(double displacement,
//     double initial_energy, double final_energy, const Vector3D& position,
//     const Vector3D& direction, const std::array<double, 4>& rnd)
// {
//     return utility_def->scattering->Scatter(
//         displacement, initial_energy, final_energy, position, direction, rnd);
// }
//
// std::pair<double, double> DirectionDeflect(CrossSection& crosssection, double particle_energy, double loss_energy)
// {
//     throw std::logic_error("Not implemented jet!");
//     /* return crosssection.StochasticDeflection(particle_energy, loss_energy); */
// }
//
// double Utility::LengthContinuous( double initial_energy, double final_energy)
// {
//     return utility_def->displacement_calc->SolveTrackIntegral(
//         initial_energy, final_energy);
// }
