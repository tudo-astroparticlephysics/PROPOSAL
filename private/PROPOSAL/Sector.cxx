
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Secondaries.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/decay/DecayChannel.h"

#include "PROPOSAL/medium/density_distr/density_distr.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/math/RandomGenerator.h"
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
          particle_,
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
          particle_,
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
      scattering_(sector.scattering_->clone(particle_, utility_))
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

    Secondaries secondaries(averaged_produced_particle_ * 1.5);

    /* std::vector<Particle*> products; */

    std::pair<double, DynamicData::Type> energy_loss;

    DynamicData continuous_loss {DynamicData::Type::None};

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

            final_energy = displacement_calculator_->GetUpperLimit(
                initial_energy, displacement_aequivaltent);
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

            final_energy = displacement_calculator_->GetUpperLimit(
                initial_energy, displacement_aequivaltent);
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
                    particle_.GetDirection(), energy_loss.second, particle_.GetEnergy(),
                    particle_.GetTime(), particle_.GetPropagatedDistance());

            /* if (energy_loss.second == DynamicData::None) { */
            /*     // in this case, no cross section is chosen, so there is no */
            /*     // interaction due to the parameterization of the cross section */
            /*     // cutoffs */
            /*     log_debug( */
            /*         "no interaction due to the parameterization of the cross " */
            /*         "section cutoffs. final energy: %f\n", */
            /*         final_energy); */
            /*     initial_energy = final_energy; */
            /*     continue; */
            /* } */

            /* if(!products.empty()){ */
            /*     // add produced particles to SecondaryVector */
            /*     for(unsigned int i=0; i<products.size(); i++){ */
            /*         products[i]->SetPosition(particle_.GetPosition()); */
            /*         products[i]->SetTime(particle_.GetTime()); */
            /*         products[i]->SetParentParticleEnergy(particle_.GetEnergy()); */
            /*     } */


            /*     if (not (sector_def_.only_loss_inside_detector && sector_def_.location != Sector::ParticleLocation::InsideDetector) ) */
            /*     { */
            /*         for (auto p : products) { */
            /*             secondaries.push_back(*p); */
            /*         } */
            /*     } */
            /* } */

            /* if(energy_loss.second != DynamicData::Particle){ */
            /*     // DynamicData::Particle means no DynamicData object will be added to SecondaryVector */
            /*     // add energy loss with DynamicData to SecondaryVector */

            /*     if (not (sector_def_.only_loss_inside_detector && sector_def_.location != Sector::ParticleLocation::InsideDetector)) */
            /*     { */
            /*             secondaries.push_back(particle_, energy_loss.second, energy_loss.first); */
            /*     } */
            /* } */

            /* if(std::get<2>(aux).second == true){ */
            /*     // fatal loss -> initial particle is destroyed */
            /*     is_decayed = true; */
            /*     final_energy -= energy_loss.first; */
            /*     break; */
            /* } */

            final_energy -= energy_loss.first;

        }else
        {
            /* products = particle_.GetDecayTable().SelectChannel().Decay(particle_); */
            /* if (not (sector_def_.only_loss_inside_detector && sector_def_.location == Sector::ParticleLocation::InsideDetector)) */
            /* { */
            /*         for (auto p : products) { */
            /*             secondaries.push_back(*p); */
            /*         } */
            /* } */

            is_decayed = true;
            final_energy = particle_.GetMass();

            // log_debug("Sampled decay of particle: %s",
            // particle_->GetName().c_str()); secondary_id    =
            // particle_->GetParticleId()  +   1;
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
        particle_.SetEnergy(particle_.GetMass());

        /* products = particle_.GetDecayTable().SelectChannel().Decay(particle_); */

        /* if (not (sector_def_.only_loss_inside_detector && sector_def_.location != Sector::ParticleLocation::InsideDetector)) */
        /* { */
        /*     for (auto p : products) { */
        /*         secondaries.push_back(*p); */
        /*     } */
        /* } */

        final_energy = particle_.GetMass();
    }

    particle_.SetEnergy(final_energy);

    averaged_produced_particle_ += (static_cast<int>(secondaries.GetNumberOfParticles())- averaged_produced_particle_ ) / n_th_call_;
    n_th_call_ += 1;

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

    scattering_->Scatter(dr, ei, ef);

    particle_.SetPropagatedDistance(dist);
    particle_.SetTime(time);
}

std::pair<double, DynamicData::Type> Sector::MakeStochasticLoss(double particle_energy)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    double total_rate = 0;
    double total_rate_weighted = 0;
    double rates_sum = 0;

    std::vector<CrossSection*> cross_sections = utility_.GetCrosssections();

    // return 0 and unknown, if there is no interaction
    std::pair<double, DynamicData::Type> energy_loss;
    energy_loss.first = 0.;
    energy_loss.second = DynamicData::None;

    /* std::pair<std::vector<Particle*>, bool> products; */

    std::vector<double> rates;

    rates.resize(cross_sections.size());

    if (sector_def_.do_stochastic_loss_weighting) {
        if (sector_def_.stochastic_loss_weighting > 0) {
            rnd2 =
                1 - rnd2 * std::pow(
                               rnd2,
                               std::abs(sector_def_.stochastic_loss_weighting));
        } else {
            rnd2 =
                rnd2 *
                std::pow(rnd2, std::abs(sector_def_.stochastic_loss_weighting));
        }
    }

    for (unsigned int i = 0; i < cross_sections.size(); i++) {
        rates[i] = cross_sections[i]->CalculatedNdx(particle_energy, rnd2);
        total_rate += rates[i];
    }

    total_rate_weighted = total_rate * rnd1;

    log_debug("Total rate = %f, total rate weighted = %f", total_rate,
              total_rate_weighted);

    for (unsigned int i = 0; i < rates.size(); i++) {
        rates_sum += rates[i];

        if (rates_sum >= total_rate_weighted) {
            energy_loss.first = cross_sections[i]->CalculateStochasticLoss(
                particle_energy, rnd2, rnd3);
            energy_loss.second = cross_sections[i]->GetTypeId();


            // Deflect the particle if necessary
            cross_sections[i]->StochasticDeflection(&particle_, particle_energy, energy_loss.first);

            // calculate produced particles, get an empty list if no particles are produced in interaction
            /* products           = cross_sections[i]->CalculateProducedParticles(particle_energy, energy_loss.first, particle_.GetDirection()); */

            break;
        }
    }


    return std::pair<double, DynamicData::Type>(energy_loss.first, energy_loss.second);
}
