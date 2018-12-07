
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/decay/DecayChannel.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/sector/Sector.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

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
    , scattering_model(ScatteringFactory::Moliere)
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
    else if (do_continuous_randomization != sector_def.do_continuous_randomization)
        return false;
    else if (do_continuous_energy_loss_output != sector_def.do_continuous_energy_loss_output)
        return false;
    else if (do_exact_time_calculation != sector_def.do_exact_time_calculation)
        return false;
    else if (scattering_model != sector_def.scattering_model)
        return false;
    else if (location != sector_def.location)
        return false;
    else if(utility_def != sector_def.utility_def)
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

Sector::Sector(Particle& particle, const Definition& sector_def)
    : sector_def_(sector_def)
    , particle_(particle)
    , geometry_(sector_def.GetGeometry().clone())
    , utility_(particle_.GetParticleDef(), sector_def.GetMedium(), sector_def.cut_settings, sector_def.utility_def)
    , displacement_calculator_(new UtilityIntegralDisplacement(utility_))
    , interaction_calculator_(new UtilityIntegralInteraction(utility_))
    , decay_calculator_(new UtilityIntegralDecay(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(ScatteringFactory::Get().CreateScattering(sector_def_.scattering_model, particle_, utility_))
{
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation)
    {
        exact_time_calculator_ = new UtilityIntegralTime(utility_);
    }

    if (sector_def_.do_continuous_randomization)
    {
        cont_rand_ = new ContinuousRandomizer(utility_);
    }
}

Sector::Sector(Particle& particle, const Definition& sector_def, const InterpolationDef& interpolation_def)
    : sector_def_(sector_def)
    , particle_(particle)
    , geometry_(sector_def.GetGeometry().clone())
    , utility_(particle_.GetParticleDef(),
               sector_def.GetMedium(),
               sector_def.cut_settings,
               sector_def.utility_def,
               interpolation_def)
    , displacement_calculator_(new UtilityInterpolantDisplacement(utility_, interpolation_def))
    , interaction_calculator_(new UtilityInterpolantInteraction(utility_, interpolation_def))
    , decay_calculator_(new UtilityInterpolantDecay(utility_, interpolation_def))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(ScatteringFactory::Get().CreateScattering(sector_def_.scattering_model,
                                                            particle_,
                                                            utility_,
                                                            interpolation_def))
{
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation)
    {
        exact_time_calculator_ = new UtilityInterpolantTime(utility_, interpolation_def);
    }

    if (sector_def_.do_continuous_randomization)
    {
        cont_rand_ = new ContinuousRandomizer(utility_, interpolation_def);
    }
}

Sector::Sector(Particle& particle, const Sector& sector)
    : sector_def_(sector.sector_def_)
    , particle_(particle)
    , geometry_(sector.geometry_->clone())
    , utility_(sector.utility_)
    , displacement_calculator_(sector.displacement_calculator_->clone(utility_))
    , interaction_calculator_(sector.interaction_calculator_->clone(utility_))
    , decay_calculator_(sector.decay_calculator_->clone(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(sector.scattering_->clone(particle_, utility_))
{
    if (particle.GetParticleDef() != sector.GetParticle().GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the sector paricle definition!");
    }

    // These are optional, therfore check NULL
    if (sector.exact_time_calculator_ != NULL)
    {
        exact_time_calculator_ = sector.exact_time_calculator_->clone(utility_);
    }

    if (sector.cont_rand_ != NULL)
    {
        cont_rand_ = new ContinuousRandomizer(utility_, *sector.cont_rand_);
    }
}

Sector::Sector(const Sector& collection)
    : sector_def_(collection.sector_def_)
    , particle_(collection.particle_)
    , geometry_(collection.geometry_->clone())
    , utility_(collection.utility_)
    , displacement_calculator_(collection.displacement_calculator_->clone(utility_))
    , interaction_calculator_(collection.interaction_calculator_->clone(utility_))
    , decay_calculator_(collection.decay_calculator_->clone(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(collection.scattering_->clone())
{
    // These are optional, therfore check NULL
    if (collection.exact_time_calculator_ != NULL)
    {
        exact_time_calculator_ = collection.exact_time_calculator_->clone(utility_);
    }

    if (collection.cont_rand_ != NULL)
    {
        cont_rand_ = new ContinuousRandomizer(utility_, *collection.cont_rand_);
    }
}

bool Sector::operator==(const Sector& sector) const
{
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
    if (exact_time_calculator_)
    {
        delete exact_time_calculator_;
    }

    if (cont_rand_)
    {
        delete cont_rand_;
    }
}

// ------------------------------------------------------------------------- //
double Sector::Propagate(double distance)
{
    bool flag;
    double displacement;

    double propagated_distance = 0;

    double initial_energy = particle_.GetEnergy();
    double final_energy   = particle_.GetEnergy();

    bool is_decayed           = false;
    bool particle_interaction = false;

    std::vector<Particle*> decay_products;

    std::pair<double, DynamicData::Type> energy_loss;

    DynamicData* continuous_loss = NULL;

    // TODO(mario): check Fri 2017/08/25
    // int secondary_id    =   0;

    // first: final energy befor first interaction second: energy at which the
    // particle decay
    // first and second are compared to decide if interaction happens or decay
    std::pair<double, double> energy_till_stochastic_;

    if (distance < 0)
    {
        distance = 0;
    }

    if (initial_energy <= particle_.GetLow() || distance == 0)
    {
        flag = false;
    } else
    {
        flag = true;
    }

    while (flag)
    {
        energy_till_stochastic_ = CalculateEnergyTillStochastic(initial_energy);

        if (energy_till_stochastic_.first > energy_till_stochastic_.second)
        {
            particle_interaction = true;
            final_energy         = energy_till_stochastic_.first;
        } else
        {
            particle_interaction = false;
            final_energy         = energy_till_stochastic_.second;
        }

        // Calculate the displacement according to initial energy and final_energy
        displacement = displacement_calculator_->Calculate(initial_energy,
                                                           final_energy,
                                                           utility_.GetMedium().GetDensityCorrection() *
                                                               (distance - propagated_distance)) /
                       utility_.GetMedium().GetDensityCorrection();

        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if (displacement > distance - propagated_distance)
        {
            displacement = distance - propagated_distance;

            final_energy = displacement_calculator_->GetUpperLimit(
                initial_energy, utility_.GetMedium().GetDensityCorrection() * displacement);
        }

        if(sector_def_.do_continuous_energy_loss_output)
        {
            continuous_loss = new DynamicData(DynamicData::ContinuousEnergyLoss);
            continuous_loss->SetPosition(particle_.GetPosition());
            continuous_loss->SetTime(particle_.GetTime());
            continuous_loss->SetParentParticleEnergy(particle_.GetEnergy());
            continuous_loss->SetPropagatedDistance(particle_.GetPropagatedDistance());
        }

        // Advance the Particle according to the displacement
        // Initial energy and final energy are needed if Molier Scattering is enabled
        AdvanceParticle(displacement, initial_energy, final_energy);

        propagated_distance += displacement;

        if (std::abs(distance - propagated_distance) < std::abs(distance) * COMPUTER_PRECISION)
        {
            propagated_distance = distance; // computer precision control
        }

        if (cont_rand_)
        {
            if (final_energy != particle_.GetLow())
            {
                final_energy =
                    cont_rand_->Randomize(initial_energy, final_energy, RandomGenerator::Get().RandomDouble());
            }
        }

        if(sector_def_.do_continuous_energy_loss_output)
        {
            continuous_loss->SetEnergy(initial_energy - final_energy);
            continuous_loss->SetDirection((particle_.GetPosition() - continuous_loss->GetPosition()));
            continuous_loss->SetTime(particle_.GetTime() - continuous_loss->GetTime());
            Output::getInstance().FillSecondaryVector(continuous_loss);
        }

        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if (final_energy == particle_.GetLow() || propagated_distance == distance)
        {
            break;
        }

        // Set the particle energy to the current energy before making
        // stochatic losses or decay
        particle_.SetEnergy(final_energy);

        if (particle_interaction)
        {
            energy_loss = MakeStochasticLoss();

            if (energy_loss.second == DynamicData::None)
            {
                // in this case, no cross section is chosen, so there is no interaction
                // due to the parameterization of the cross section cutoffs
                log_debug("no interaction due to the parameterization of the cross section cutoffs. final energy: %f\n",
                          final_energy);
                initial_energy = final_energy;
                continue;
            }
            final_energy -= energy_loss.first;
            Output::getInstance().FillSecondaryVector(particle_, energy_loss.second, energy_loss.first);

            // log_debug("Energyloss: %f\t%s", energy_loss.first,
            // secondary_id    =   particle_.GetParticleId() + 1;
        } else
        {
            decay_products = particle_.GetDecayTable().SelectChannel().Decay(particle_);
            Output::getInstance().FillSecondaryVector(decay_products);

            is_decayed   = true;
            final_energy = particle_.GetMass();

            // log_debug("Sampled decay of particle: %s", particle_->GetName().c_str());
            // secondary_id    = particle_->GetParticleId()  +   1;
        }

        // break if the lower limit of particle energy is reached
        if (final_energy <= particle_.GetLow())
        {
            break;
        }

        // Next round: update the inital energy
        initial_energy = final_energy;
    }

    // if a particle is below a specific energy 'elow' and stopping_decay is enabled,
    // the muon gets forced to decay, instead of propagating all the time with dEdx
    // with no significantly produced light
    if (sector_def_.stopping_decay && propagated_distance != distance && !is_decayed)
    {
        // TODO: understand what happens in the two following lines
        //       why is the particle energy set to its mass?
        //       why is the time increased randomly?
        //       it doesn't make sense for me (jsoedingrekso)
        // particle_.GetEnergy() = particle_.GetParticleDef().mass;
        // particle_.GetTime() -= particle_.GetLifetime()*std::log(RandomGenerator::Get().RandomDouble());
        decay_products = particle_.GetDecayTable().SelectChannel().Decay(particle_);
        Output::getInstance().FillSecondaryVector(decay_products);

        final_energy = particle_.GetMass();
    }

    particle_.SetEnergy(final_energy);

    // Particle reached the border, final energy is returned
    if (propagated_distance == distance)
    {
        return final_energy;
    }
    // The particle stopped/decayed, the propageted distance is return with a minus sign
    else
    {
        return -propagated_distance;
    }
}

std::pair<double, double> Sector::CalculateEnergyTillStochastic(double initial_energy)
{
    double rndd = -std::log(RandomGenerator::Get().RandomDouble());
    double rndi = -std::log(RandomGenerator::Get().RandomDouble());

    double rndiMin = 0;
    double rnddMin = 0;

    std::pair<double, double> final;

    // solving the tracking integral
    if (particle_.GetLifetime() < 0)
    {
        rnddMin = 0;
    } else
    {
        rnddMin = decay_calculator_->Calculate(initial_energy, particle_.GetParticleDef().low, rndd) /
                  utility_.GetMedium().GetDensityCorrection();
    }

    rndiMin = interaction_calculator_->Calculate(initial_energy, particle_.GetParticleDef().low, rndi);

    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0)
    {
        final.second = particle_.GetLow();
    } else
    {
        final.second =
            decay_calculator_->GetUpperLimit(initial_energy, rndd * utility_.GetMedium().GetDensityCorrection());
    }

    if (rndi >= rndiMin || rndiMin <= 0)
    {
        final.first = particle_.GetLow();
    } else
    {
        final.first = interaction_calculator_->GetUpperLimit(initial_energy, rndi);
    }

    return final;
}

void Sector::AdvanceParticle(double dr, double ei, double ef)
{

    double dist       = particle_.GetPropagatedDistance();
    double time       = particle_.GetTime();
    Vector3D position = particle_.GetPosition();

    dist += dr;

    if (exact_time_calculator_)
    {
        time += exact_time_calculator_->Calculate(ei, ef, 0.0) / utility_.GetMedium().GetDensityCorrection();
    } else
    {
        time += dr / SPEED;
    }

    scattering_->Scatter(dr, ei, ef);

    particle_.SetPropagatedDistance(dist);
    particle_.SetTime(time);
}

std::pair<double, DynamicData::Type> Sector::MakeStochasticLoss()
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    double total_rate          = 0;
    double total_rate_weighted = 0;
    double rates_sum           = 0;

    std::vector<CrossSection*> cross_sections = utility_.GetCrosssections();

    // return 0 and unknown, if there is no interaction
    std::pair<double, DynamicData::Type> energy_loss;
    energy_loss.first  = 0.;
    energy_loss.second = DynamicData::None;

    std::vector<double> rates;

    rates.resize(cross_sections.size());

    if (sector_def_.do_stochastic_loss_weighting)
    {
        if (sector_def_.stochastic_loss_weighting > 0)
        {
            rnd2 = 1 - rnd2 * std::pow(rnd2, std::abs(sector_def_.stochastic_loss_weighting));
        }
        else
        {
            rnd2 = rnd2 * std::pow(rnd2, std::abs(sector_def_.stochastic_loss_weighting));
        }
    }

    for (unsigned int i = 0; i < cross_sections.size(); i++)
    {
        rates[i] = cross_sections[i]->CalculatedNdx(particle_.GetEnergy(), rnd2);
        total_rate += rates[i];
    }

    total_rate_weighted = total_rate * rnd1;

    log_debug("Total rate = %f, total rate weighted = %f", total_rate, total_rate_weighted);

    for (unsigned int i = 0; i < rates.size(); i++)
    {
        rates_sum += rates[i];

        if (rates_sum > total_rate_weighted)
        {
            energy_loss.first  = cross_sections[i]->CalculateStochasticLoss(particle_.GetEnergy(), rnd2, rnd3);
            energy_loss.second = cross_sections[i]->GetTypeId();
            break;
        }
    }

    return energy_loss;
}
