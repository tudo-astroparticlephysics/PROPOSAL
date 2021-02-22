
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/Spherical3D.h"

using namespace PROPOSAL;

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

PropagationUtility::PropagationUtility(
    PropagationUtility::Collection const& collect)
    : collection(collect)
{
    if (collect.interaction_calc == nullptr
        or collect.displacement_calc == nullptr
        or collect.time_calc == nullptr) {
        throw std::invalid_argument("Interaction, displacement and time "
                                    "calculator need to be defined.");
    }
}

Interaction::Loss PropagationUtility::EnergyStochasticloss(double energy,
                                                           double rnd)
{
    auto rates = collection.interaction_calc->Rates(energy);
    auto loss = collection.interaction_calc->SampleLoss(energy, rates, rnd);

    return loss;
}

double PropagationUtility::EnergyDecay(
    double energy, std::function<double()> rnd, double density)
{
    if (collection.decay_calc) {
        return collection.decay_calc->EnergyDecay(energy, rnd(), density);
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
    double initial_energy, double final_energy, double distance, double density)
{
    return collection.time_calc->TimeElapsed(
        initial_energy, final_energy, distance, density);
}

std::tuple<Cartesian3D, Cartesian3D> PropagationUtility::DirectionsScatter(
    double displacement, double initial_energy, double final_energy,
    const Vector3D& direction, std::function<double()> rnd)
{
    if (collection.scattering) {
        std::array<double, 4> random_numbers;
        for (auto& r : random_numbers)
            r = rnd();
        auto random_angles = collection.scattering->CalculateMultipleScattering(
            displacement, initial_energy, final_energy, random_numbers);

        auto& sx = random_angles[multiple_scattering::Parametrization::SX];
        auto& sy = random_angles[multiple_scattering::Parametrization::SY];
        auto& tx = random_angles[multiple_scattering::Parametrization::TX];
        auto& ty = random_angles[multiple_scattering::Parametrization::TY];
        auto sz = std::sqrt(std::max(1. - (sx * sx + sy * sy), 0.));
        auto tz = std::sqrt(std::max(1. - (tx * tx + ty * ty), 0.));

        auto direction_spherical = Spherical3D(direction);
        auto sinth = std::sin(direction_spherical.GetZenith());
        auto costh = std::cos(direction_spherical.GetZenith());
        auto sinph = std::sin(direction_spherical.GetAzimuth());
        auto cosph = std::cos(direction_spherical.GetAzimuth());

        auto rotate_vector_x = Cartesian3D(costh * cosph, costh * sinph, -sinth);
        auto rotate_vector_y = Cartesian3D(-sinph, cosph, 0.);

        // Rotation towards all tree axes
        auto direction_cartesian = Cartesian3D(direction);
        auto mean_direction = sz * direction_cartesian;
        mean_direction += sx * rotate_vector_x;
        mean_direction += sy * rotate_vector_y;

        // Rotation towards all tree axes
        auto final_direction = tz * direction_cartesian;
        final_direction += tx * rotate_vector_x;
        final_direction += ty * rotate_vector_y;

        return std::make_tuple(mean_direction, final_direction);
    }
    auto dir = Cartesian3D(direction.GetCartesianCoordinates());
    return std::make_tuple(dir, dir); // no scattering
}

Cartesian3D PropagationUtility::DirectionDeflect(InteractionType type,
    double initial_energy, double final_energy, const Vector3D& direction,
    std::function<double()> rnd) const
{
    if (collection.scattering) {
        auto v_rnd = std::vector<double>(
            collection.scattering->StochasticDeflectionRandomNumbers(type));
        for (auto& r : v_rnd)
            r = rnd();
        auto angles = collection.scattering->CalculateStochasticDeflection(
            type, initial_energy, final_energy, v_rnd);
        auto direction_new = Cartesian3D(direction);
        direction_new.deflect(std::cos(angles[0]), angles[1]);
    }
    return direction;
}

double PropagationUtility::LengthContinuous(
    double initial_energy, double final_energy)
{
    return collection.displacement_calc->SolveTrackIntegral(
        initial_energy, final_energy);
}
