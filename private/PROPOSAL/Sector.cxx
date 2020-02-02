
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Secondaries.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/decay/DecayChannel.h"

#include "PROPOSAL/medium/density_distr/density_distr.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Sector.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include <memory>
using namespace PROPOSAL;

/******************************************************************************
 *                                 Sector                                 *
 ******************************************************************************/

Sector::Definition::Definition()
    : do_stochastic_loss_weighting(false),
      stochastic_loss_weighting(0),
      stopping_decay(true),
      do_continuous_randomization(true),
      do_continuous_energy_loss_output(false),
      do_exact_time_calculation(true),
      only_loss_inside_detector(false),
      scattering_model(ScatteringFactory::HighlandIntegral),
      location(Sector::ParticleLocation::InsideDetector),
      utility_def(),
      cut_settings(),
      medium_(new Ice()),
      geometry_(new Sphere(Vector3D(), 1.0e20, 0.0)) {}

Sector::Definition::Definition(const Definition& def)
    : do_stochastic_loss_weighting(def.do_stochastic_loss_weighting),
      stochastic_loss_weighting(def.stochastic_loss_weighting),
      stopping_decay(def.stopping_decay),
      do_continuous_randomization(def.do_continuous_randomization),
      do_continuous_energy_loss_output(def.do_continuous_energy_loss_output),
      do_exact_time_calculation(def.do_exact_time_calculation),
      only_loss_inside_detector(def.only_loss_inside_detector),
      scattering_model(def.scattering_model),
      location(def.location),
      utility_def(def.utility_def),
      cut_settings(def.cut_settings),
      medium_(def.medium_->clone()),
      geometry_(def.geometry_->clone()) {}

bool Sector::Definition::operator==(const Definition& sector_def) const {
    if (do_stochastic_loss_weighting != sector_def.do_stochastic_loss_weighting)
        return false;
    else if (stochastic_loss_weighting != sector_def.stochastic_loss_weighting)
        return false;
    else if (stopping_decay != sector_def.stopping_decay)
        return false;
    else if (do_continuous_randomization !=
             sector_def.do_continuous_randomization)
        return false;
    else if (do_continuous_energy_loss_output !=
             sector_def.do_continuous_energy_loss_output)
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

bool Sector::Definition::operator!=(const Definition& sector_def) const {
    return !(*this == sector_def);
}

void Sector::Definition::swap(Definition& definition) {
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
Sector::Definition& Sector::Definition::operator=(
    const Definition& definition) {
    if (this != &definition) {
        Sector::Definition tmp(definition);
        swap(tmp);
    }
    return *this;
}

Sector::Definition::~Definition() {
    delete medium_;
    delete geometry_;
}

void Sector::Definition::SetMedium(const Medium& medium) {
    delete medium_;
    medium_ = medium.clone();
}

void Sector::Definition::SetGeometry(const Geometry& geometry) {
    delete geometry_;
    geometry_ = geometry.clone();
}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

Sector::Sector(Particle& particle, const Definition& sector_def)
    : sector_def_(sector_def),
      particle_(particle),
      geometry_(sector_def.GetGeometry().clone()),
      utility_(particle_.GetParticleDef(),
               sector_def.GetMedium(),
               sector_def.cut_settings,
               sector_def.utility_def),
      displacement_calculator_(new UtilityIntegralDisplacement(utility_)),
      interaction_calculator_(new UtilityIntegralInteraction(utility_)),
      decay_calculator_(new UtilityIntegralDecay(utility_)),
      exact_time_calculator_(NULL),
      cont_rand_(NULL),
      scattering_(ScatteringFactory::Get().CreateScattering(
          sector_def_.scattering_model,
          particle_.GetParticleDef(),
          utility_)) {
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation) {
        exact_time_calculator_ = new UtilityIntegralTime(utility_);
    }

    if (sector_def_.do_continuous_randomization) {
        cont_rand_ = new ContinuousRandomizer(utility_);
    }
}

Sector::Sector(Particle& particle,
               const Definition& sector_def,
               const InterpolationDef& interpolation_def)
    : sector_def_(sector_def),
      particle_(particle),
      geometry_(sector_def.GetGeometry().clone()),
      utility_(particle_.GetParticleDef(),
               sector_def.GetMedium(),
               sector_def.cut_settings,
               sector_def.utility_def,
               interpolation_def),
      displacement_calculator_(
          new UtilityInterpolantDisplacement(utility_, interpolation_def)),
      interaction_calculator_(
          new UtilityInterpolantInteraction(utility_, interpolation_def)),
      decay_calculator_(
          new UtilityInterpolantDecay(utility_, interpolation_def)),
      exact_time_calculator_(NULL),
      cont_rand_(NULL),
      scattering_(ScatteringFactory::Get().CreateScattering(
          sector_def_.scattering_model,
          particle_.GetParticleDef(),
          utility_,
          interpolation_def)) {
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation) {
        exact_time_calculator_ =
            new UtilityInterpolantTime(utility_, interpolation_def);
    }

    if (sector_def_.do_continuous_randomization) {
        cont_rand_ = new ContinuousRandomizer(utility_, interpolation_def);
    }
}

Sector::Sector(Particle& particle, const Sector& sector)
    : sector_def_(sector.sector_def_),
      particle_(particle),
      geometry_(sector.geometry_->clone()),
      utility_(sector.utility_),
      displacement_calculator_(sector.displacement_calculator_->clone(utility_)),
      interaction_calculator_(sector.interaction_calculator_->clone(utility_)),
      decay_calculator_(sector.decay_calculator_->clone(utility_)),
      exact_time_calculator_(NULL),
      cont_rand_(NULL),
      scattering_(sector.scattering_->clone(particle_.GetParticleDef(), utility_))
{
    if (particle.GetParticleDef() != sector.GetParticle().GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the sector particle definition!");
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
    : sector_def_(sector.sector_def_),
      particle_(sector.particle_),
      geometry_(sector.geometry_->clone()),
      utility_(sector.utility_),
      displacement_calculator_(
          sector.displacement_calculator_->clone(utility_)),
      interaction_calculator_(sector.interaction_calculator_->clone(utility_)),
      decay_calculator_(sector.decay_calculator_->clone(utility_)),
      exact_time_calculator_(NULL),
      cont_rand_(NULL),
      scattering_(sector.scattering_->clone()) {
    // These are optional, therfore check NULL
    if (sector.exact_time_calculator_ != NULL) {
        exact_time_calculator_ = sector.exact_time_calculator_->clone(utility_);
    }

    if (sector.cont_rand_ != NULL) {
        cont_rand_ = new ContinuousRandomizer(utility_, *sector.cont_rand_);
    }
}

bool Sector::operator==(const Sector& sector) const {
    if (sector_def_ != sector.sector_def_)
        return false;
    else if (particle_ != sector.particle_)
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

bool Sector::operator!=(const Sector& sector) const {
    return !(*this == sector);
}

Sector::~Sector() {
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
std::pair<double, Secondaries> Sector::Propagate(double distance) {

    bool flag;
    double displacement;

    double propagated_distance = 0;

    double initial_energy = particle_.GetEnergy();
    double final_energy = particle_.GetEnergy();

    bool is_decayed = false;
    bool particle_interaction = false;

    Secondaries secondaries;
    secondaries.reserve(static_cast<size_t>(produced_particle_moments_.first + 2 * std::sqrt(produced_particle_moments_.second)));

    /* std::vector<Particle*> products; */

    std::pair<double, int> energy_loss;

    DynamicData continuous_loss {0};

    // TODO(mario): check Fri 2017/08/25
    // int secondary_id    =   0;

    // first: final energy before first interaction second: energy at which the
    // particle decay
    // first and second are compared to decide if interaction happens or decay
    std::pair<double, double> energy_till_stochastic_;

    if (distance < 0) {
        distance = 0;
    }

    if (initial_energy <= particle_.GetLow() || distance == 0) {
        flag = false;
    } else {
        flag = true;
    }

    while (flag) {
        energy_till_stochastic_ = CalculateEnergyTillStochastic(initial_energy);

        if (energy_till_stochastic_.first > energy_till_stochastic_.second) {
            particle_interaction = true;
            final_energy = energy_till_stochastic_.first;
        } else {
            particle_interaction = false;
            final_energy = energy_till_stochastic_.second;
        }

        try {
            displacement = displacement_calculator_->Calculate(
                initial_energy, final_energy, distance - propagated_distance,
                particle_.GetPosition(), particle_.GetDirection());
        } catch (DensityException& e) {
            displacement = distance - propagated_distance;

            double displacement_aequivaltent =
                utility_.GetMedium().GetDensityDistribution().Calculate(
                    particle_.GetPosition(), particle_.GetDirection(),
                    displacement);

            final_energy = EnergyDistance(initial_energy, displacement_aequivaltent);
        }

        // The first interaction or decay happens behind the distance we want to
        // propagate So we calculate the final energy using only continuous
        // losses
        if (displacement > distance - propagated_distance) {
            displacement = distance - propagated_distance;

            double displacement_aequivaltent =
                utility_.GetMedium().GetDensityDistribution().Calculate(
                    particle_.GetPosition(), particle_.GetDirection(),
                    displacement);
            final_energy = EnergyDistance(initial_energy, displacement_aequivaltent);
        }

        if (sector_def_.do_continuous_energy_loss_output) {
            continuous_loss.SetPosition(particle_.GetPosition());
            continuous_loss.SetTime(particle_.GetTime());
            continuous_loss.SetParentParticleEnergy(particle_.GetEnergy());
            continuous_loss.SetPropagatedDistance(
                particle_.GetPropagatedDistance());
        }

        // Advance the Particle according to the displacement
        // Initial energy and final energy are needed if Molier Scattering is
        // enabled
        AdvanceParticle(displacement, initial_energy, final_energy);

        propagated_distance += displacement;

        if (std::abs(distance - propagated_distance) <
            std::abs(distance) * COMPUTER_PRECISION) {
            propagated_distance = distance;  // computer precision control
        }

        if (cont_rand_) {
            if (final_energy != particle_.GetLow()) {
                final_energy = cont_rand_->Randomize(
                    initial_energy, final_energy,
                    RandomGenerator::Get().RandomDouble());
            }
        }

        if (sector_def_.do_continuous_energy_loss_output) {
            continuous_loss.SetEnergy(initial_energy - final_energy);
            continuous_loss.SetDirection(particle_.GetDirection());
            continuous_loss.SetPropagatedDistance( particle_.GetPropagatedDistance() - continuous_loss.GetPropagatedDistance());
            if (not (sector_def_.only_loss_inside_detector && sector_def_.location != Sector::ParticleLocation::InsideDetector)) {
                    secondaries.push_back(continuous_loss);
            }
        }

        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if (final_energy == particle_.GetLow() ||
            propagated_distance == distance) {
            break;
        }

        // Set the particle energy to the current energy before making
        // stochastic losses or decay
        particle_.SetEnergy(final_energy);

        if (particle_interaction){
            energy_loss = MakeStochasticLoss(final_energy);

            secondaries.emplace_back(energy_loss.second, particle_.GetPosition(),
                    particle_.GetDirection(), energy_loss.first, particle_.GetEnergy(),
                    particle_.GetTime(), particle_.GetPropagatedDistance());

            final_energy -= energy_loss.first;

        }else
        {
            is_decayed = true;
            final_energy = particle_.GetMass();

            secondaries.emplace_back(static_cast<int>(InteractionType::Decay), particle_.GetPosition(),
                    particle_.GetDirection(), final_energy, particle_.GetEnergy(),
                    particle_.GetTime(), particle_.GetPropagatedDistance());
        }

        // break if the lower limit of particle energy is reached
        if (final_energy <= particle_.GetLow()) {
            break;
        }

        // Next round: update the initial energy
        initial_energy = final_energy;
    }

    // if a particle is below a specific energy 'elow' and stopping_decay is
    // enabled, the muon gets forced to decay, instead of propagating all the
    // time with dEdx with no significantly produced light
    if (sector_def_.stopping_decay && propagated_distance != distance &&
        !is_decayed) {
        // The time is shifted due to the exponential lifetime.
        double particle_time = particle_.GetTime();
        particle_time -= particle_.GetLifetime() *
                         std::log(RandomGenerator::Get().RandomDouble());
        particle_.SetTime(particle_time);

        // TODO: one should also advance the particle according to the sampeled
        // time and set the new position as the endpoint.

        final_energy = particle_.GetMass();

        secondaries.emplace_back(static_cast<int>(InteractionType::Decay), particle_.GetPosition(),
                particle_.GetDirection(), final_energy, particle_.GetEnergy(),
                particle_.GetTime(), particle_.GetPropagatedDistance());

        particle_.SetEnergy(particle_.GetMass());
    }

    particle_.SetEnergy(final_energy);

    n_th_call_ += 1.;
    double produced_particles_ = static_cast<double>(secondaries.GetNumberOfParticles());
    produced_particle_moments_ = welfords_online_algorithm(produced_particles_, n_th_call_, produced_particle_moments_.first, produced_particle_moments_.second);

    // Particle reached the border, final energy is returned
    if (propagated_distance == distance) {
        return std::make_pair(final_energy, secondaries);
    }
    // The particle stopped/decayed, the propagated distance is return with a
    // minus sign
    else {
        return std::make_pair(-propagated_distance, secondaries);
    }
}

std::pair<double, double> Sector::CalculateEnergyTillStochastic(
    double initial_energy) {
    double rndd = -std::log(RandomGenerator::Get().RandomDouble());
    double rndi = -std::log(RandomGenerator::Get().RandomDouble());

    double rndiMin = 0;
    double rnddMin = 0;

    std::pair<double, double> final;

    // solving the tracking integral
    if (particle_.GetLifetime() < 0) {
        rnddMin = 0;
    } else {
        rnddMin = decay_calculator_->Calculate(
                      initial_energy, particle_.GetParticleDef().low, rndd) /
                  utility_.GetMedium().GetDensityDistribution().Evaluate(
                      particle_.GetPosition());
    }

    rndiMin = interaction_calculator_->Calculate(
        initial_energy, particle_.GetParticleDef().low, rndi);

    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0) {
        final.second = particle_.GetLow();
    } else {
        final.second = decay_calculator_->GetUpperLimit(
            initial_energy,
            rndd * utility_.GetMedium().GetDensityDistribution().Evaluate(
                       particle_.GetPosition()));
    }

    if (rndi >= rndiMin || rndiMin <= 0) {
        final.first = particle_.GetLow();
    } else {
        final.first =
            interaction_calculator_->GetUpperLimit(initial_energy, rndi);
    }

    return final;
}

double Sector::EnergyDistance(double initial_energy, double disp) {
    return displacement_calculator_->GetUpperLimit(initial_energy, disp);
};

void Sector::AdvanceParticle(double dr, double ei, double ef) {
    double dist = particle_.GetPropagatedDistance();
    double time = particle_.GetTime();

    dist += dr;

    if (exact_time_calculator_) {
        // DensityDistribution Approximation: Use the DensityDistribution at the
        // position of initial energy
        time += exact_time_calculator_->Calculate(ei, ef, 0.0) /
                utility_.GetMedium().GetDensityDistribution().Evaluate(
                    particle_.GetPosition());
    } else {
        time += dr / SPEED;
    }


    if( sector_def_.scattering_model != ScatteringFactory::Enum::NoScattering ){
        Directions directions = scattering_->Scatter(dr, ei, ef, particle_.GetPosition(), particle_.GetDirection());
        particle_.SetPosition(particle_.GetPosition() + dr * directions.u_);
        particle_.SetDirection(directions.n_i_);
    } else {
        particle_.SetPosition(particle_.GetPosition() + dr * particle_.GetDirection());
    }


    particle_.SetPropagatedDistance(dist);
    particle_.SetTime(time);
}

std::pair<double, int> Sector::MakeStochasticLoss(double particle_energy)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    // do an optional bias to the randomly sampled energy loss
    if (sector_def_.do_stochastic_loss_weighting) {
        double aux = rnd2 * std::pow(rnd2, std::abs(sector_def_.stochastic_loss_weighting));
        if (sector_def_.stochastic_loss_weighting > 0) {
            rnd2 = 1 - aux;
        } else {
            rnd2 = aux;
        }
    }

    std::pair<double, int> energy_loss;
    std::pair<std::vector<Particle*>, bool> products;
    std::pair<double, double> deflection_angles;

    energy_loss = utility_.StochasticLoss(particle_energy, rnd1, rnd2, rnd3);

    std::vector<CrossSection*> cross_sections = utility_.GetCrosssections();
    for (unsigned int i = 0; i < cross_sections.size(); i++) {
        if (cross_sections[i]->GetTypeId() == energy_loss.second)
        {
            deflection_angles = cross_sections[i]->StochasticDeflection(particle_energy, energy_loss.first);
            // TODO: dirty hack until alternative output branch is merged and particle deflection is calculated later
            if (deflection_angles.first != particle_energy)
            {
                particle_.DeflectDirection(deflection_angles.first, deflection_angles.second);
            }
            // calculate produced particles, get an empty list if no particles are produced in interaction
            products = cross_sections[i]->CalculateProducedParticles(particle_energy, energy_loss.first, particle_.GetDirection());
            break;
        }

    }

    return std::pair<double, int>(energy_loss.first, energy_loss.second);
}
