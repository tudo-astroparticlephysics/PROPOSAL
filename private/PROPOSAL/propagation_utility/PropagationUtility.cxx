
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;
using std::tuple;

/*
PropagationUtility::Definition::Definition(CrossSectionList cross,
    const ParticleDef& p_def, std::shared_ptr<Scattering> scattering = nullptr,
    std::shared_ptr<InterpolationDef> inter_def = nullptr)
    : scattering(scattering)
{
    if (!inter_def) {
        log_warn("No interpolation definition defined. Integral will not be "
                 "approximate by interpolants. Performance will be poor.");

        displacement_calc.reset(new
DisplacementBuilder<UtilityIntegral>(cross)); interaction_calc.reset(new
InteractionBuilder<UtilityIntegral>(cross)); decay_calc.reset(new
DecayBuilder<UtilityIntegral>(cross, p_def)); } else {
        displacement_calc.reset(new
DisplacementBuilder<UtilityInterpolant>(cross)); interaction_calc.reset(new
InteractionBuilder<UtilityInterpolant>(cross)); decay_calc.reset(new
DecayBuilder<UtilityInterpolant>(cross, p_def));
    }

    if (!scattering) {
        log_debug("No scattering defined. Particle will only be deflected in "
                  "stochastic interactions");
    }

    if (!cont_rand) {
        log_debug("No continuous randomization used.");
    } else {
        if (!inter_def) {
            cont_rand.reset(new ContRandBuilder<UtilityIntegral>(cross));
        } else {
            cont_rand.reset(new ContRandBuilder<UtilityInterpolant>(cross));
        }
    }

    if (!time_calc) {
        log_debug("Exact time calculation disabled: Velocity of particles will "
                  "be approximated as speed of light");
    } else {
        if (!inter_def) {
            time_calc.reset(new ExactTimeBuilder<UtilityIntegral>(cross,
p_def)); } else { time_calc.reset(new
ExactTimeBuilder<UtilityInterpolant>(cross, p_def));
        }
    }
}

*/

bool PropagationUtility::Collection::operator==(const Collection& lhs)
{
    if (interaction_calc != lhs.interaction_calc)
        return false;
    if (displacement_calc != lhs.displacement_calc)
        return false;
    if (time_calc != lhs.time_calc)
        return false;
    if (scattering != lhs.scattering)
        return false;
    if (decay_calc != lhs.decay_calc)
        return false;
    if (cont_rand != lhs.cont_rand)
        return false;
    return true;
}

PropagationUtility::PropagationUtility(PropagationUtility::Collection collect)
    : collection(collect)
{
    if (collect.interaction_calc == nullptr
        or collect.displacement_calc == nullptr
        or collect.time_calc == nullptr) {
        throw std::invalid_argument("Interaction, displacement and time "
                                    "calculator need to be defined.");
    }
}

/*
PropagationUtility::PropagationUtility(const PropagationUtility& utility) :
utility_def(utility.utility_def){
}
*/

tuple<InteractionType, const Component*, double> PropagationUtility::EnergyStochasticloss(
    double energy, double rnd)
{
    return collection.interaction_calc->TypeInteraction(energy, rnd);
}

double PropagationUtility::EnergyDecay(
    double energy, std::function<double()> rnd)
{
    if (collection.decay_calc) {
        return collection.decay_calc->EnergyDecay(energy, rnd());
    }
    return 0; // no decay, e.g. particle is stable
}

double PropagationUtility::EnergyInteraction(
    double energy, std::function<double()> rnd)
{
    return collection.interaction_calc->EnergyInteraction(energy, rnd());
}

double PropagationUtility::EnergyRandomize(
    double initial_energy, double final_energy, std::function<double()> rnd)
{
    if (collection.cont_rand) {
        final_energy = collection.cont_rand->EnergyRandomize(
            initial_energy, final_energy, rnd());
    }
    return final_energy; // no randomization
}

double PropagationUtility::EnergyDistance(
    double initial_energy, double distance)
{
    return collection.displacement_calc->UpperLimitTrackIntegral(
        initial_energy, distance);
}

double PropagationUtility::TimeElapsed(
    double initial_energy, double final_energy, double distance)
{
    return collection.time_calc->TimeElapsed(
        initial_energy, final_energy, distance);
}

tuple<Vector3D, Vector3D> PropagationUtility::DirectionsScatter(
    double displacement, double initial_energy, double final_energy,
    const Vector3D& direction, const std::array<double, 4>& rnd)
{
    if (collection.scattering) {
        return collection.scattering->Scatter(
            displacement, initial_energy, final_energy, direction, rnd);
    }
    return std::make_tuple(direction, direction); // no scattering
}

/* std::pair<double, double> DirectionDeflect( */
/*     CrossSection& crosssection, double particle_energy, double loss_energy)
 */
/* { */
/*     return crosssection.StochasticDeflection(particle_energy, loss_energy);
 */
/* } */

double PropagationUtility::LengthContinuous(
    double initial_energy, double final_energy)
{
    return collection.displacement_calc->SolveTrackIntegral(
        initial_energy, final_energy);
}
