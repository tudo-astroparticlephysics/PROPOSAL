
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Secondaries.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/decay/DecayChannel.h"

#include "PROPOSAL/medium/density_distr/density_distr.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/Sector.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include <memory>
#include <tuple>
using namespace PROPOSAL;

/******************************************************************************
 *                                 Sector                                 *
 ******************************************************************************/

Sector::Definition::Definition()
    : do_stochastic_loss_weighting(false)
    , stochastic_loss_weighting(0)
    , stopping_decay(true)
    , do_continuous_randomization(true)
    , do_continuous_energy_loss_output(false)
    , do_exact_time_calculation(true)
    , only_loss_inside_detector(false)
    , scattering_model(ScatteringFactory::HighlandIntegral)
    , location(Sector::ParticleLocation::InsideDetector)
    , utility_def()
    , cut_settings()
    , medium_(new Ice())
    , geometry_(new Sphere(Vector3D(), 1.0e20, 0.0))
{
}

Sector::Definition::Definition(const Definition& def)
    : do_stochastic_loss_weighting(def.do_stochastic_loss_weighting)
    , stochastic_loss_weighting(def.stochastic_loss_weighting)
    , stopping_decay(def.stopping_decay)
    , do_continuous_randomization(def.do_continuous_randomization)
    , do_continuous_energy_loss_output(def.do_continuous_energy_loss_output)
    , do_exact_time_calculation(def.do_exact_time_calculation)
    , only_loss_inside_detector(def.only_loss_inside_detector)
    , scattering_model(def.scattering_model)
    , location(def.location)
    , utility_def(def.utility_def)
    , cut_settings(def.cut_settings)
    , medium_(def.medium_->clone())
    , geometry_(def.geometry_->clone())
{
}

bool Sector::Definition::operator==(const Definition& sector_def) const
{
    if (do_stochastic_loss_weighting != sector_def.do_stochastic_loss_weighting)
        return false;
    else if (stochastic_loss_weighting != sector_def.stochastic_loss_weighting)
        return false;
    else if (stopping_decay != sector_def.stopping_decay)
        return false;
    else if (do_continuous_randomization
        != sector_def.do_continuous_randomization)
        return false;
    else if (do_continuous_energy_loss_output
        != sector_def.do_continuous_energy_loss_output)
        return false;
    else if (do_exact_time_calculation != sector_def.do_exact_time_calculation)
        return false;
    else if (only_loss_inside_detector != sector_def.only_loss_inside_detector)
        return false;
    else if (scattering_model != sector_def.scattering_model)
        return false;
    else if (location != sector_def.location)
        return false;
    else if (utility_def != sector_def.utility_def)
        return false;
    else if (*medium_ != *sector_def.medium_)
        return false;
    else if (*geometry_ != *sector_def.geometry_)
        return false;
    return true;
}

bool Sector::Definition::operator!=(const Definition& sector_def) const
{
    return !(*this == sector_def);
}

void Sector::Definition::swap(Definition& definition)
{
    using std::swap;

    swap(do_stochastic_loss_weighting, definition.do_stochastic_loss_weighting);
    swap(stochastic_loss_weighting, definition.stochastic_loss_weighting);
    swap(stopping_decay, definition.stopping_decay);
    swap(do_continuous_randomization, definition.do_continuous_randomization);
    swap(do_continuous_energy_loss_output,
        definition.do_continuous_energy_loss_output);
    swap(do_exact_time_calculation, definition.do_exact_time_calculation);
    swap(only_loss_inside_detector, definition.only_loss_inside_detector);
    swap(scattering_model, definition.scattering_model);
    swap(location, definition.location);
    swap(utility_def, definition.utility_def);
    medium_->swap(*definition.medium_);
    geometry_->swap(*definition.geometry_);
}
Sector::Definition& Sector::Definition::operator=(const Definition& definition)
{
    if (this != &definition) {
        Sector::Definition tmp(definition);
        swap(tmp);
    }
    return *this;
}

Sector::Definition::~Definition()
{
    delete medium_;
    delete geometry_;
}

void Sector::Definition::SetMedium(const Medium& medium)
{
    delete medium_;
    medium_ = medium.clone();
}

void Sector::Definition::SetGeometry(const Geometry& geometry)
{
    delete geometry_;
    geometry_ = geometry.clone();
}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

Sector::Sector(const ParticleDef& particle_def, const Definition& sector_def)
    : sector_def_(sector_def)
    , particle_def_(particle_def)
    , geometry_(sector_def.GetGeometry().clone())
    , utility_(particle_def, sector_def.GetMedium(), sector_def.cut_settings,
          sector_def.utility_def)
    , displacement_calculator_(new UtilityIntegralDisplacement(utility_))
    , interaction_calculator_(new UtilityIntegralInteraction(utility_))
    , decay_calculator_(new UtilityIntegralDecay(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(ScatteringFactory::Get().CreateScattering(
          sector_def_.scattering_model, particle_def, utility_))
{
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation) {
        exact_time_calculator_ = new UtilityIntegralTime(utility_);
    }

    if (sector_def_.do_continuous_randomization) {
        cont_rand_ = new ContinuousRandomizer(utility_);
    }
}

Sector::Sector(const ParticleDef& particle_def, const Definition& sector_def,
    const InterpolationDef& interpolation_def)
    : sector_def_(sector_def)
    , particle_def_(particle_def)
    , geometry_(sector_def.GetGeometry().clone())
    , utility_(particle_def, sector_def.GetMedium(), sector_def.cut_settings,
          sector_def.utility_def, interpolation_def)
    , displacement_calculator_(
          new UtilityInterpolantDisplacement(utility_, interpolation_def))
    , interaction_calculator_(
          new UtilityInterpolantInteraction(utility_, interpolation_def))
    , decay_calculator_(
          new UtilityInterpolantDecay(utility_, interpolation_def))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(ScatteringFactory::Get().CreateScattering(
          sector_def_.scattering_model, particle_def, utility_,
          interpolation_def))
{
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation) {
        exact_time_calculator_
            = new UtilityInterpolantTime(utility_, interpolation_def);
    }

    if (sector_def_.do_continuous_randomization) {
        cont_rand_ = new ContinuousRandomizer(utility_, interpolation_def);
    }
}

Sector::Sector(const ParticleDef& particle_def, const Sector& sector)
    : sector_def_(sector.sector_def_)
    , particle_def_(particle_def)
    , geometry_(sector.geometry_->clone())
    , utility_(sector.utility_)
    , displacement_calculator_(sector.displacement_calculator_->clone(utility_))
    , interaction_calculator_(sector.interaction_calculator_->clone(utility_))
    , decay_calculator_(sector.decay_calculator_->clone(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(sector.scattering_->clone(particle_def, utility_))
{
    if (particle_def_ != sector.GetParticleDef()) {
        log_fatal("Particle definition should be equal to the sector particle "
                  "definition!");
    }

    // These are optional, therfore check NULL
    if (sector.exact_time_calculator_ != NULL) {
        exact_time_calculator_ = sector.exact_time_calculator_->clone(utility_);
    }

    if (sector.cont_rand_ != NULL) {
        cont_rand_ = new ContinuousRandomizer(utility_, *sector.cont_rand_);
    }
}

Sector::Sector(const Sector& sector)
    : sector_def_(sector.sector_def_)
    , particle_def_(sector.particle_def_)
    , geometry_(sector.geometry_->clone())
    , utility_(sector.utility_)
    , displacement_calculator_(sector.displacement_calculator_->clone(utility_))
    , interaction_calculator_(sector.interaction_calculator_->clone(utility_))
    , decay_calculator_(sector.decay_calculator_->clone(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(sector.scattering_->clone())
{
    // These are optional, therfore check NULL
    if (sector.exact_time_calculator_ != NULL) {
        exact_time_calculator_ = sector.exact_time_calculator_->clone(utility_);
    }

    if (sector.cont_rand_ != NULL) {
        cont_rand_ = new ContinuousRandomizer(utility_, *sector.cont_rand_);
    }
}

bool Sector::operator==(const Sector& sector) const
{
    if (sector_def_ != sector.sector_def_)
        return false;
    else if (particle_def_ != sector.particle_def_)
        return false;
    else if (*geometry_ != *sector.geometry_)
        return false;
    else if (utility_ != sector.utility_)
        return false;
    else if (*cont_rand_ != *sector.cont_rand_)
        return false;
    else if (*scattering_ != *sector.scattering_)
        return false;
    return true;
}

bool Sector::operator!=(const Sector& sector) const
{
    return !(*this == sector);
}

Sector::~Sector()
{
    delete geometry_;
    delete scattering_;

    delete displacement_calculator_;
    delete interaction_calculator_;
    delete decay_calculator_;

    // These are optional, therfore check NULL
    if (exact_time_calculator_) {
        delete exact_time_calculator_;
    }

    if (cont_rand_) {
        delete cont_rand_;
    }
}

// ------------------------------------------------------------------------- //
double Sector::DecayEnergy(double initial_energy, double rnd)
{
    double rndd = -std::log(rnd);
    double rnddMin = 0;

    std::pair<double, double> final;

    // solving the tracking integral
    if (particle_def_.lifetime < 0) {
        return particle_def_.low;
    }

    rnddMin
        = decay_calculator_->Calculate(initial_energy, particle_def_.low, rndd);

    // evaluating the energy loss
    if (rndd < rnddMin || rnddMin <= 0) {
        return particle_def_.low;
    }

    return decay_calculator_->GetUpperLimit(initial_energy, rndd);
}

double Sector::InteractionEnergy(double initial_energy, double rnd)
{
    double rndi = -std::log(rnd);
    double rndiMin = 0;

    // solving the tracking integral
    rndiMin = interaction_calculator_->Calculate(
        initial_energy, particle_def_.low, rndi);

    if (rndi >= rndiMin || rndiMin <= 0) {
        return particle_def_.low;
    }

    return interaction_calculator_->GetUpperLimit(initial_energy, rndi);
}

double Sector::BorderEnergy(double distance, double rnd)
{
    return displacement_calculator_->GetUpperLimit(distance, rnd);
}

double Sector::GetDisplacement(const DynamicData p_condition,
    const double final_energy, const double border_length)
{
    return displacement_calculator_->Calculate(p_condition.GetEnergy(),
        final_energy, border_length, p_condition.GetPosition(),
        p_condition.GetDirection());
}

double Sector::BorderLength(const DynamicData p_condition)
{
    // loop ueber alle sektoren hoeherer ordnung die im aktuellen sektor liegen
    // und getroffen werden koennten
    return geometry_
        ->DistanceToBorder(
            p_condition.GetPosition(), p_condition.GetDirection())
        .first;
}

std::vector<std::tuple<int, double, double>> Sector::GetSteplengths(
    const DynamicData p_condition)
{
    double rnd1{ RandomGenerator::Get().RandomDouble() };
    double rnd2{ RandomGenerator::Get().RandomDouble() };

    double distance_border{ BorderLength(p_condition) };

    double energy_interaction{ InteractionEnergy(
        p_condition.GetEnergy(), rnd1) };
    double energy_decay{ DecayEnergy(p_condition.GetEnergy(), rnd2) };

    double distance_interaction{ GetDisplacement(
        p_condition, energy_interaction, distance_border) };
    double distance_decay{ GetDisplacement(
        p_condition, energy_decay, distance_border) };

    int InteractionTypeId = 0;
    int DecayTypeId = 1;
    int BorderTypeId = 2;

    double borderEnergy = 0;

    auto interaction = std::make_tuple(
        InteractionTypeId, energy_interaction, distance_interaction);
    auto decay = std::make_tuple(DecayTypeId, energy_decay, distance_decay);
    auto border = std::make_tuple(BorderTypeId, borderEnergy, distance_border);

    std::vector<std::tuple<int, double, double>> stepLengths{ interaction,
        decay, border };

    return stepLengths;
}

std::tuple<int, double, double> Sector::minimizeSteplengths(
    std::vector<std::tuple<int, double, double>> stepLengths)
{
    std::tuple<int, double, double> minimalStepLength = stepLengths.front();
    for (auto& i : stepLengths) {
        if (std::get<2>(i) <= std::get<2>(minimalStepLength)) {
            minimalStepLength = i;
        }
    }
    return minimalStepLength;
}

std::shared_ptr<DynamicData> Sector::DoInteraction(
    const DynamicData p_condition)
{
    std::pair<int, double> stochastic_loss
        = MakeStochasticLoss(p_condition.GetEnergy());

    CrossSection* cross_section
        = utility_.GetCrosssection(stochastic_loss.second);
    std::pair<double, double> deflection_angles
        = cross_section->StochasticDeflection(
            p_condition.GetEnergy(), stochastic_loss.first);

    Vector3D new_direction(p_condition.GetDirection());
    new_direction.deflect(deflection_angles.first, deflection_angles.second);

    double new_energy = p_condition.GetEnergy() - stochastic_loss.first;
    return std::make_shared<DynamicData>(stochastic_loss.second,
        p_condition.GetPosition(), new_direction, new_energy,
        p_condition.GetEnergy(), p_condition.GetTime(),
        p_condition.GetPropagatedDistance());
}

std::shared_ptr<DynamicData> Sector::DoDecay(const DynamicData p_condition)
{

    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::Decay), p_condition.GetPosition(),
        p_condition.GetDirection(), p_condition.GetEnergy(),
        p_condition.GetParentParticleEnergy(), p_condition.GetTime(),
        p_condition.GetPropagatedDistance());
}

double Sector::CalculateTime(
    const DynamicData& p_condition, double final_energy, double displacement)
{
    if (exact_time_calculator_) {
        // DensityDistribution Approximation: Use the DensityDistribution at the
        // position of initial energy
        return p_condition.GetTime()
            + exact_time_calculator_->Calculate(
                  p_condition.GetEnergy(), final_energy, 0.0)
            / utility_.GetMedium().GetDensityDistribution().Evaluate(
                  p_condition.GetPosition());
    }

    return p_condition.GetTime() + displacement / SPEED;
}

std::shared_ptr<DynamicData> Sector::DoBorder(const DynamicData p_condition)
{
    double distance = BorderLength(p_condition);
    double final_energy = BorderEnergy(distance, 0);
    double time = CalculateTime(p_condition, final_energy, distance);
    Vector3D position
        = p_condition.GetPosition() + distance * p_condition.GetDirection();

    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::ContinuousEnergyLoss), position,
        p_condition.GetDirection(), final_energy, p_condition.GetEnergy(), time,
        distance);
}

std::shared_ptr<DynamicData> Sector::DoContinuous(
    const DynamicData p_condition, double final_energy, double displacement)
{
    double dist{ p_condition.GetPropagatedDistance() };
    double time{ p_condition.GetTime() };

    Vector3D position{ p_condition.GetPosition() };
    Vector3D direction{ p_condition.GetDirection() };

    dist += displacement;

    time = CalculateTime(p_condition, final_energy, displacement);

    if (sector_def_.scattering_model != ScatteringFactory::Enum::NoScattering) {
        Directions directions = scattering_->Scatter(displacement,
            p_condition.GetEnergy(), final_energy, position, direction);
        position = position + displacement * directions.u_;
        direction = directions.n_i_;
    } else {
        position = position + displacement * direction;
    }

    if (cont_rand_) {
        if (final_energy != particle_def_.low) {
            final_energy = cont_rand_->Randomize(p_condition.GetEnergy(),
                final_energy, RandomGenerator::Get().RandomDouble());
        }
    }

    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::ContinuousEnergyLoss), position,
        direction, final_energy, p_condition.GetEnergy(), time, dist);
}

std::pair<double, int> Sector::MakeStochasticLoss(double particle_energy)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    std::pair<double, int> energy_loss;
    std::pair<double, double> deflection_angles;

    // do an optional bias to the randomly sampled energy loss
    if (sector_def_.do_stochastic_loss_weighting) {
        double aux = rnd2
            * std::pow(rnd2, std::abs(sector_def_.stochastic_loss_weighting));
        if (sector_def_.stochastic_loss_weighting > 0) {
            rnd2 = 1 - aux;
        } else {
            rnd2 = aux;
        }
    }

    energy_loss = utility_.StochasticLoss(particle_energy, rnd1, rnd2, rnd3);

    return energy_loss;
}

Secondaries Sector::Propagate(
    DynamicData& p_initial, double distance, double minimal_energy)
{

    distance += p_initial.GetPropagatedDistance();
    std::shared_ptr<DynamicData> p_condition
        = std::make_shared<DynamicData>(p_initial);
    Secondaries secondaries;

    do {
        secondaries.push_back(*p_condition);

        std::vector<std::tuple<int, double, double>> steplengths
            = GetSteplengths(*p_condition);
        std::tuple<int, double, double> steplength
            = minimizeSteplengths(steplengths);

        auto continuous = DoContinuous(
            *p_condition, std::get<1>(steplength), std::get<2>(steplength));
        if (std::get<0>(steplength) == 0) {
            p_condition = std::move(DoInteraction(*continuous));
        }
        if (std::get<0>(steplength) == 1) {
            p_condition = std::move(DoDecay(*continuous));
            secondaries.push_back(*p_condition);
            break;
        }
        if (std::get<0>(steplength) == 2) {
            p_condition = std::move(DoBorder(*continuous));
            secondaries.push_back(*p_condition);
            break;
        }

        std::cout << "p_condition->GetPropagatedDistance():"
                  << p_condition->GetPropagatedDistance() << std::endl;
        std::cout << "p_condition->GetEnergy():" << p_condition->GetEnergy()
                  << std::endl;
        std::cout << "distance > p_condition->GetPropagatedDistance():"
                  << (distance > p_condition->GetPropagatedDistance())
                  << std::endl;
        std::cout << "minimal_energy < p_condition->GetEnergy()"
                  << (minimal_energy < p_condition->GetEnergy()) << std::endl;
    } while (distance > p_condition->GetPropagatedDistance()
        and minimal_energy < p_condition->GetEnergy());

    DynamicData p_last = secondaries.GetSecondaries().back();

    double border_distance{ p_last.GetPropagatedDistance() };
    double border_energy{ p_last.GetEnergy() };

    DynamicData p_final(
        static_cast<int>(InteractionType::ContinuousEnergyLoss));
    if (distance < border_distance) {
        border_distance = distance - p_last.GetPropagatedDistance();
        border_energy = BorderEnergy(border_distance, border_distance);
    }

    if (minimal_energy > border_energy) {
        border_energy = minimal_energy;
        border_distance
            = GetDisplacement(p_last, minimal_energy, border_distance);
        distance =p_last.GetPropagatedDistance() + border_distance;
    }
    p_final.SetPosition(
        p_last.GetPosition() + border_distance * p_last.GetDirection());
    p_final.SetDirection(p_last.GetDirection());

    p_final.SetParentParticleEnergy(p_last.GetEnergy());
    p_final.SetTime(CalculateTime(p_last, p_last.GetEnergy(), border_distance));
    p_final.SetPropagatedDistance(distance);
    p_final.SetEnergy(border_energy);

    secondaries.push_back(p_final);

    return secondaries;
}
