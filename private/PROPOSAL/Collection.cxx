#include <boost/bind.hpp>

#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Collection.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Decay.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Geometry.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/ScatteringDefault.h"
#include "PROPOSAL/ScatteringMoliere.h"
#include "PROPOSAL/ScatteringFirstOrder.h"
#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                               CollectionDef                                *
******************************************************************************/

CollectionDef::CollectionDef()
    : do_weighting(false)
    , weighting_order(0)
    , do_scattering(false)
    , do_continuous_randomization_(false)
    , lpm_effect_enabled(false)
    , do_exact_time_calculation(false)
    , location(0)
    , density_correction(1.0)
    , order_of_interpolation(5)
    , raw(true)
    , path_to_tables("")
{
}

CollectionDef::~CollectionDef()
{
}

/******************************************************************************
*                                 Collection                                 *
******************************************************************************/


// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

// Standard constructor
Collection::Collection()
    : ini_(0)
    , collection_def_()
    , weighting_starts_at_(0)
    , geometry_(new Sphere(Vector3D(), 1e18, 0))
    , medium_(new Water())
    , cut_settings_()
    , randomizer_(NULL)
    , scattering_(NULL)
    , crosssections_()
{
    crosssections_.push_back(new Bremsstrahlung(medium_, &cut_settings_));
    crosssections_.push_back(new Epairproduction(medium_, &cut_settings_));
    crosssections_.push_back(new Photonuclear(medium_, &cut_settings_));
    crosssections_.push_back(new Ionization(medium_, &cut_settings_));

    //TODO(mario): Polymorphic initilaization in collections childs  Sun 2017/08/27
    if (collection_def_.do_continuous_randomization_)
    {
        randomizer_ = new ContinuousRandomization();
    }

    //TODO(mario): Polymorphic initilaization in collections childs  Sun 2017/08/27
    if (collection_def_.do_scattering)
    {
        scattering_ = new ScatteringDefault();
    }
}

Collection::Collection(const Medium& medium,
                       const Geometry& geometry,
                       const EnergyCutSettings& cut_settings,
                       const CollectionDef& def)
    : ini_(0)
    , collection_def_(def)
    , weighting_starts_at_(0)
    , geometry_(geometry.clone())
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , randomizer_(NULL)
    , scattering_(NULL)
    , crosssections_()
{
    crosssections_.push_back(new Bremsstrahlung(medium_, &cut_settings_));
    crosssections_.push_back(new Ionization(medium_, &cut_settings_));
    crosssections_.push_back(new Photonuclear(medium_, &cut_settings_));
    crosssections_.push_back(new Epairproduction(medium_, &cut_settings_));

    //TODO(mario): Polymorphic initilaization in collections childs  Sun 2017/08/27
    if (collection_def_.do_continuous_randomization_)
    {
        randomizer_ = new ContinuousRandomization();
    }

    //TODO(mario): Polymorphic initilaization in collections childs  Sun 2017/08/27
    if (collection_def_.do_scattering)
    {
        scattering_ = new ScatteringDefault();
    }
}

Collection::Collection(const Collection& collection)
    :ini_(collection.ini_)
    ,collection_def_(collection.collection_def_)
    ,weighting_starts_at_(collection.weighting_starts_at_)
    ,geometry_(collection.geometry_->clone())
    ,medium_(collection.medium_->clone())
    ,cut_settings_(collection.cut_settings_)
    ,randomizer_(collection.randomizer_) //TODO(mario): ranomizer clone Sat 2017/08/26
    ,scattering_(collection.scattering_) //TODO(mario): scatter clone Sat 2017/08/26
{
    crosssections_.resize(collection.crosssections_.size());

    //TODO(mario): clone Sat 2017/08/26
    for(unsigned int i =0; i<collection.crosssections_.size(); i++)
    {
        switch (collection.crosssections_.at(i)->GetType())
        {
            case ParticleType::Brems:
                crosssections_.at(i) = new Bremsstrahlung( *(Bremsstrahlung*)collection.crosssections_.at(i) );
                break;
            case ParticleType::DeltaE:
                crosssections_.at(i) = new Ionization( *(Ionization*)collection.crosssections_.at(i) );
                break;
            case ParticleType::EPair:
                crosssections_.at(i) = new Epairproduction( *(Epairproduction*)collection.crosssections_.at(i) );
                break;
            case ParticleType::NuclInt:
                crosssections_.at(i) = new Photonuclear( *(Photonuclear*)collection.crosssections_.at(i) );
                break;
            default:
                log_fatal("Unknown cross section");
                exit(1);
        }
    }
}

Collection::~Collection()
{
    delete medium_;
    delete geometry_;

    if (randomizer_)
    {
        delete randomizer_;
    }

    //TODO(mario): delete scatter Sat 2017/08/26
}

// ------------------------------------------------------------------------- //
double Collection::Propagate(PROPOSALParticle& particle, double distance)
{
    bool flag;
    double displacement;

    double propagated_distance = 0;

    double initial_energy = particle.GetEnergy();
    double final_energy   = particle.GetEnergy();

    bool particle_interaction = false;

    pair<double, ParticleType::Enum> decay;
    std::vector<PROPOSALParticle*> decay_products;

    pair<double, ParticleType::Enum> energy_loss;

    // TODO(mario): check Fri 2017/08/25
    // int secondary_id    =   0;

    // first: final energy befor first interaction second: energy at which the
    // particle decay
    // first and second are compared to decide if interaction happens or decay
    pair<double, double> energy_till_stochastic_;

    if (distance < 0)
    {
        distance = 0;
    }

    if (initial_energy <= particle.GetLow() || distance == 0)
    {
        flag = false;
    } else
    {
        flag = true;
    }

    while (flag)
    {
        energy_till_stochastic_ = CalculateEnergyTillStochastic(particle, initial_energy);
        if (energy_till_stochastic_.first > energy_till_stochastic_.second)
        {
            particle_interaction = true;
            final_energy         = energy_till_stochastic_.first;
        } else
        {
            particle_interaction = false;
            final_energy         = energy_till_stochastic_.second;
        }

        // Calculate the displacement according to initial energy initial_energy and final_energy
        displacement =
            CalculateDisplacement(
                particle, initial_energy, final_energy, collection_def_.density_correction * (distance - propagated_distance)) /
            collection_def_.density_correction;

        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if (displacement > distance - propagated_distance)
        {
            displacement = distance - propagated_distance;

            final_energy = CalculateFinalEnergy(particle, initial_energy, collection_def_.density_correction * displacement);
        }
        // Advance the Particle according to the displacement
        // Initial energy and final energy are needed if Molier Scattering is enabled
        AdvanceParticle(particle, displacement, initial_energy, final_energy);

        propagated_distance += displacement;

        if (abs(distance - propagated_distance) < abs(distance) * COMPUTER_PRECISION)
        {
            propagated_distance = distance; // computer precision control
        }

        //TODO(mario): Revert randomizer Fri 2017/08/25
        // Randomize the continuous energy loss if this option is enabled
        if (collection_def_.do_continuous_randomization_)
        {
            if (final_energy != particle.GetLow())
            {
                double rnd = RandomGenerator::Get().RandomDouble();
                final_energy = randomizer_->Randomize(particle, crosssections_, initial_energy, final_energy, rnd);
            }
        }

        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if (final_energy == particle.GetLow() || propagated_distance == distance)
        {
            break;
        }

        // Set the particle energy to the current energy before making
        // stochatic losses or decay
        particle.SetEnergy(final_energy);

        if (particle_interaction)
        {
            energy_loss = MakeStochasticLoss(particle);
            if (energy_loss.second == ParticleType::unknown)
            {
                // in this case, no cross section is chosen, so there is no interaction
                // due to the parameterization of the cross section cutoffs
                log_debug("no interaction due to the parameterization of the cross section cutoffs. final energy: %f\n",
                          final_energy);
                initial_energy = final_energy;
                continue;
            }
            final_energy -= energy_loss.first;
            // log_debug("Energyloss: %f\t%s", energy_loss.first,
            // PROPOSALParticle::GetName(energy_loss.second).c_str());
            // //TODO(mario): hack Thu 2017/08/24
            Output::getInstance().FillSecondaryVector(&particle, ParticleDef(BremsDef::Get()), energy_loss.first, 0);
            // secondary_id    =   particle.GetParticleId() + 1;
            // Output::getInstance().FillSecondaryVector(&particle, secondary_id, energy_loss, 0);
        } else
        {
            DecayChannel* mode = &particle.GetDecayTable().SelectChannel();
            decay_products     = mode->Decay(&particle);
            Output::getInstance().FillSecondaryVector(decay_products);

            // TODO(mario): Delete decay products Tue 2017/08/22

            // decay           =   current_collection_->MakeDecay();
            // final_energy    =   0;
            // log_debug("Decay of particle: %s", particle_->GetName().c_str());
            // secondary_id    = particle_->GetParticleId()  +   1;
            // Output::getInstance().FillSecondaryVector(particle_, secondary_id, decay ,0);
        }

        // break if the lower limit of particle energy is reached
        if (final_energy <= particle.GetLow())
        {
            break;
        }

        // Next round: update the inital energy
        initial_energy = final_energy;
    }

    // if(stopping_decay_)
    // {
    //     if(propagated_distance!=distance && final_energy!=0 && particle_->GetLifetime()>=0)
    //     {
    //         particle_->SetEnergy(particle_->GetMass());
    //
    //         double t    =   particle_->GetT() -particle_->GetLifetime()*log(RandomDouble());
    //         double product_energy   =   0;
    //
    //         pair<double, ParticleType::Enum> decay_to_store;
    //         secondary_id    =   particle_->GetParticleId() + 1;
    //
    //         particle_->SetT( t );
    //
    //         if(particle_->GetType()==2)
    //         {
    //             // --------------------------------------------------------------------- //
    //             // Calculate random numbers before passing to a fuction, because
    //             // the order of argument evaluation is unspecified in c++ standards and
    //             // therfor depend on the compiler.
    //             // --------------------------------------------------------------------- //
    //
    //             double rnd1 = RandomDouble();
    //             double rnd2 = RandomDouble();
    //
    //             product_energy  =   current_collection_->GetDecay()->CalculateProductEnergy(rnd1, 0.5, rnd2);
    //         }
    //         else
    //         {
    //             product_energy  =   current_collection_->GetDecay()->CalculateProductEnergy(RandomDouble(), 0.5,
    //             0.5);
    //         }
    //
    //         decay_to_store.first    =   product_energy;
    //         decay_to_store.second   =   current_collection_->GetDecay()->GetOut();
    //
    //         final_energy  =   0;
    //
    //         Output::getInstance().FillSecondaryVector(particle_,secondary_id, decay_to_store, final_energy);
    //     }
    // }

    particle.SetEnergy(final_energy);

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
    // Should never be here
    return 0;
}

std::pair<double, double> Collection::CalculateEnergyTillStochastic(const PROPOSALParticle& particle,
                                                                    double initial_energy)
{
    double rndd = -log(RandomGenerator::Get().RandomDouble());
    double rndi = -log(RandomGenerator::Get().RandomDouble());

    double rndiMin = 0;
    double rnddMin = 0;

    pair<double, double> final;

    // solving the tracking integral
    if (particle.GetLifetime() < 0)
    {
        rnddMin = 0;
    } else
    {
        rnddMin = CalculateTrackingIntegal(particle, initial_energy, rndd, false) / collection_def_.density_correction;
    }

    rndiMin = CalculateTrackingIntegal(particle, initial_energy, rndi, true);
    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0)
    {
        final.second = particle.GetLow();
    } else
    {
        final.second = CalculateFinalEnergy(particle, initial_energy, rndd * collection_def_.density_correction, false);
    }

    if (rndi >= rndiMin || rndiMin <= 0)
    {
        final.first = particle.GetLow();
    } else
    {
        final.first = CalculateFinalEnergy(particle, initial_energy, rndi, true);
    }

    return final;
}

void Collection::AdvanceParticle(PROPOSALParticle& particle, double dr, double ei, double ef)
{

    double dist       = particle.GetPropagatedDistance();
    double time       = particle.GetT();
    Vector3D position = particle.GetPosition();

    dist += dr;

    if (collection_def_.do_exact_time_calculation)
    {
        time += CalculateParticleTime(particle, ei, ef) / collection_def_.density_correction;
    } else
    {
        time += dr / SPEED;
    }

    // TODO(mario): Adjucst the whole scatteing class Thu 2017/08/24
    if (collection_def_.do_scattering)
    {
        scattering_->Scatter(particle, crosssections_, dr, ei, ef);
    }

    // if(scattering_model_!=-1)
    // {
    //     switch(scattering_model_)
    //     {
    //         case 0:
    //             current_collection_->GetScattering()->Scatter(dr,ei,ef);
    //             break;
    //
    //         case 1:
    //             scatteringFirstOrder_->Scatter(dr, particle_, current_collection_->GetMedium());
    //             break;
    //
    //         case 2:
    //             scatteringFirstOrderMoliere_->Scatter(dr, particle_, current_collection_->GetMedium());
    //             break;
    //         default:
    //             log_error("Never should be here! scattering_model = %i !",scattering_model_);
    //     }
    //
    // }
    // else
    // {
    //     position = position + dr*particle.GetDirection();
    //     particle.SetPosition(position);
    // }

    particle.SetPropagatedDistance(dist);
    particle.SetT(time);
}

pair<double, ParticleType::Enum> Collection::MakeStochasticLoss(const PROPOSALParticle& particle)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    double total_rate          = 0;
    double total_rate_weighted = 0;
    double rates_sum           = 0;

    // return 0 and unknown, if there is no interaction
    pair<double, ParticleType::Enum> energy_loss;
    energy_loss.first  = 0.;
    energy_loss.second = ParticleType::unknown;

    std::vector<double> rates;

    rates.resize(crosssections_.size());

    if (collection_def_.do_weighting)
    {
        if (particle.GetPropagatedDistance() > weighting_starts_at_)
        {
            double exp   = abs(collection_def_.weighting_order);
            double power = pow(rnd2, exp);

            if (collection_def_.weighting_order > 0)
            {
                rnd2 = 1 - power * rnd2;
            } else
            {
                rnd2 = power * rnd2;
            }

            collection_def_.weighting_order     = (1 + exp) * power;
            weighting_starts_at_ = particle.GetPropagatedDistance();
            collection_def_.do_weighting        = false;
        }
    }
    // if (particle_->GetEnergy() < 650) printf("energy: %f\n", particle_->GetEnergy());
    for (unsigned int i = 0; i < GetCrosssections().size(); i++)
    {
        rates.at(i) = crosssections_.at(i)->CalculatedNdx(particle, rnd2);
        total_rate += rates.at(i);
        // if (rates.at(i) == 0) printf("%i = 0, energy: %f\n", i, particle_->GetEnergy());
        log_debug("Rate for %s = %f", crosssections_.at(i)->GetName().c_str(), rates.at(i));
    }

    total_rate_weighted = total_rate * rnd1;

    log_debug("Total rate = %f, total rate weighted = %f", total_rate, total_rate_weighted);

    for (unsigned int i = 0; i < rates.size(); i++)
    {
        rates_sum += rates.at(i);

        if (rates_sum > total_rate_weighted)
        {
            energy_loss.first  = crosssections_.at(i)->CalculateStochasticLoss(particle, rnd2, rnd3);
            energy_loss.second = crosssections_.at(i)->GetType();
            break;
        }
    }

    return energy_loss;
}

double Collection::MakeDecay(const PROPOSALParticle& particle)
{
    // TODO(mario): multiplier? Was not used before Fri 2017/08/25
    // if(multiplier_ <= 0 || particle.GetLifetime() < 0)
    // {
    //     return 0;
    // }
    //
    // return multiplier_/max((particle.GetMomentum()/particle.GetMass())*particle.GetLifetime()*SPEED, XRES);
    if (particle.GetLifetime() < 0)
    {
        return 0;
    }

    return max((particle.GetMomentum() / particle.GetMass()) * particle.GetLifetime() * SPEED, XRES);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// DecayChannel::DecayProducts Collection::MakeDecay()
// {
//     const DecayChannel* mode = particle_->GetDecayTable().SelectChannel();
//     return mode->Decay(particle_);
// }

// pair<double, ParticleType::Enum> Collection::MakeDecay()
// {
//     // --------------------------------------------------------------------- //
//     // Calculate random numbers before passing to a fuction, because
//     // the order of argument evaluation is unspecified in c++ standards and
//     // therfor depend on the compiler.
//     // --------------------------------------------------------------------- //
//
//     pair<double, ParticleType::Enum> decay_pair;
//     const DecayChannel* mode = particle_->GetDecayTable().SelectChannel();
//     decay_pair.first = mode->Decay(particle_);
//
//     // double rnd1 = MathModel::RandomDouble();
//     // double rnd2 = MathModel::RandomDouble();
//     // double rnd3 = MathModel::RandomDouble();
//     //
//     // return MakeDecay(rnd1, rnd2, rnd3);
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// pair<double, ParticleType::Enum> Collection::MakeDecay(double rnd1,double rnd2, double rnd3)
// {
//     pair<double, ParticleType::Enum> decay;
//
//     if(particle_->GetType() == ParticleType::TauPlus || particle_->GetType() == ParticleType::TauMinus)
//     {
//         decay.first     =   decay_->CalculateProductEnergy(rnd1, rnd2, rnd3);
//     }
//     else
//     {
//         decay.first     =   decay_->CalculateProductEnergy(rnd2, rnd3, 0.5);
//     }
//
//     decay.second    =   decay_->GetOut();
//
//     return decay;
// }

// double Collection::Randomize(double initial_energy, double final_energy)
// {
//     double rnd = RandomGenerator::Get().RandomDouble();
//     return randomizer_->Randomize(initial_energy, final_energy, rnd);
// }

// ------------------------------------------------------------------------- //
// Lpm effect & randomization
// ------------------------------------------------------------------------- //


void Collection::EnableLpmEffect()
{
    collection_def_.lpm_effect_enabled = true;
    for (unsigned int i = 0; i < crosssections_.size(); i++)
    {
        crosssections_.at(i)->EnableLpmEffect(collection_def_.lpm_effect_enabled);
    }
}

void Collection::DisableLpmEffect()
{
    collection_def_.lpm_effect_enabled = false;
    for (unsigned int i = 0; i < crosssections_.size(); i++)
    {
        crosssections_.at(i)->EnableLpmEffect(collection_def_.lpm_effect_enabled);
    }
}

// void Collection::EnableContinuousRandomization()
// {
//     randomizer_                  = new ContinuousRandomization(particle_, medium_, crosssections_);
//     do_continuous_randomization_ = true;
// }
//
// void Collection::DisableContinuousRandomization()
// {
//     delete randomizer_;
//     randomizer_                  = NULL;
//     do_continuous_randomization_ = false;
// }

void Collection::EnableExactTimeCalculation()
{
    collection_def_.do_exact_time_calculation   =   true;
}

void Collection::DisableExactTimeCalculation()
{
    collection_def_.do_exact_time_calculation   =   false;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// Copyconstructor
// Collection::Collection(const Collection &collection)
//     : MathModel(collection)
//     ,order_of_interpolation_     ( collection.order_of_interpolation_ )
//     ,do_interpolation_           ( collection.do_interpolation_ )
//     ,lpm_effect_enabled_         ( collection.lpm_effect_enabled_ )
//     ,ini_                        ( collection.ini_ )
//     ,do_weighting_               ( collection.do_weighting_ )
//     ,weighting_order_            ( collection.weighting_order_ )
//     ,weighting_starts_at_        ( collection.weighting_starts_at_ )
//     ,enable_randomization_       ( collection.enable_randomization_ )
//     ,do_continuous_randomization_( collection.do_continuous_randomization_ )
//     ,do_scattering_              ( collection.do_scattering_ )
//     ,location_                   ( collection.location_ )
//     ,density_correction_         ( collection.density_correction_ )
//     ,do_time_interpolation_      ( collection.do_time_interpolation_ )
//     ,do_exact_time_calculation_   ( collection.do_exact_time_calculation_ )
//     ,up_                         ( collection.up_)
//     ,bigLow_                     ( collection.bigLow_ )
//     ,storeDif_                   ( collection.storeDif_ )
//     ,particle_                   ( new PROPOSALParticle(*collection.particle_) )
//     ,backup_particle_            ( new PROPOSALParticle(*collection.backup_particle_) )
//     ,medium_                     ( collection.medium_->clone() )
//     ,integral_                   ( new Integral(*collection.integral_) )
//     ,cut_settings_               ( new EnergyCutSettings(*collection.cut_settings_) )
//     ,decay_                      ( new Decay(*collection.decay_) )
//     ,prop_decay_                 ( new Integral(*collection.prop_decay_) )
//     ,prop_interaction_           ( new Integral(*collection.prop_interaction_) )
//     ,time_particle_              ( new Integral(*collection.time_particle_) )
//
// {
//     crosssections_.resize(collection.crosssections_.size());
//
//     for(unsigned int i =0; i<collection.crosssections_.size(); i++)
//     {
//         switch (collection.crosssections_.at(i)->GetType())
//         {
//             case ParticleType::Brems:
//                 crosssections_.at(i) = new Bremsstrahlung( *(Bremsstrahlung*)collection.crosssections_.at(i) );
//                 break;
//             // case ParticleType::DeltaE:
//             //     crosssections_.at(i) = new Ionization( *(Ionization*)collection.crosssections_.at(i) );
//             //     break;
//             // case ParticleType::EPair:
//             //     crosssections_.at(i) = new Epairproduction( *(Epairproduction*)collection.crosssections_.at(i) );
//             //     break;
//             // case ParticleType::NuclInt:
//             //     crosssections_.at(i) = new Photonuclear( *(Photonuclear*)collection.crosssections_.at(i) );
//             //     break;
//             default:
//                 log_fatal("Unknown cross section");
//                 exit(1);
//         }
//     }
//
//     if(collection.interpolant_ != NULL)
//     {
//         interpolant_ = new Interpolant(*collection.interpolant_) ;
//     }
//     else
//     {
//         interpolant_ = NULL;
//     }
//
//     if(collection.interpolant_diff_ != NULL)
//     {
//         interpolant_diff_ = new Interpolant(*collection.interpolant_diff_) ;
//     }
//     else
//     {
//         interpolant_diff_ = NULL;
//     }
//
//     if(collection.interpol_prop_decay_ != NULL)
//     {
//         interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_) ;
//     }
//     else
//     {
//         interpol_prop_decay_ = NULL;
//     }
//
//     if(collection.interpol_prop_decay_diff_ != NULL)
//     {
//         interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_) ;
//     }
//     else
//     {
//         interpol_prop_decay_diff_ = NULL;
//     }
//
//     if(collection.interpol_prop_interaction_ != NULL)
//     {
//         interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_) ;
//     }
//     else
//     {
//         interpol_prop_interaction_ = NULL;
//     }
//
//     if(collection.interpol_prop_interaction_diff_ != NULL)
//     {
//         interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_) ;
//     }
//     else
//     {
//         interpol_prop_interaction_diff_ = NULL;
//     }
//
//     if(collection.interpol_time_particle_ != NULL)
//     {
//         interpol_time_particle_ = new Interpolant(*collection.interpol_time_particle_) ;
//     }
//     else
//     {
//         interpol_time_particle_ = NULL;
//     }
//
//     if(collection.interpol_time_particle_diff_ != NULL)
//     {
//         interpol_time_particle_diff_ = new Interpolant(*collection.interpol_time_particle_diff_) ;
//     }
//     else
//     {
//         interpol_time_particle_diff_ = NULL;
//     }
//
//     if(collection.randomizer_ != NULL)
//     {
//         randomizer_ = new ContinuousRandomization(*collection.randomizer_) ;
//     }
//     else
//     {
//         randomizer_ = NULL;
//     }
//
//     if(collection.scattering_ != NULL)
//     {
//         scattering_ = new Scattering(*collection.scattering_) ;
//     }
//     else
//     {
//         scattering_ = NULL;
//     }
//
//     if(collection.geometry_ != NULL)
//     {
//         // geometry_ = new Geometry(*collection.geometry_) ;
//         geometry_ = collection.geometry_->clone();
//     }
//     else
//     {
//         geometry_ = NULL;
//     }
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// Collection& Collection::operator=(const Collection &collection)
// {
//     if (this != &collection)
//     {
//       Collection tmp(collection);
//       swap(tmp);
//     }
//     return *this;
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// bool Collection::operator==(const Collection &collection) const
// {
//     if( order_of_interpolation_     != collection.order_of_interpolation_ )  return false;
//     if( do_interpolation_           != collection.do_interpolation_ )        return false;
//     if( lpm_effect_enabled_         != collection.lpm_effect_enabled_ )      return false;
//     if( ini_                        != collection.ini_ )                     return false;
//     if( *particle_                  != *collection.particle_ )               return false;
//     if( *medium_                    != *collection.medium_ )                 return false;
//     if( *integral_                  != *collection.integral_ )               return false;
//     if( *cut_settings_              != *collection.cut_settings_ )           return false;
//     if( *prop_decay_                != *collection.prop_decay_ )             return false;
//     if( *prop_interaction_          != *collection.prop_interaction_ )       return false;
//     if( up_                         != collection.up_)                       return false;
//     if( do_weighting_               != collection.do_weighting_ )            return false;
//     if( weighting_order_            != collection.weighting_order_ )         return false;
//     if( weighting_starts_at_        != collection.weighting_starts_at_ )     return false;
//     if( enable_randomization_       != collection.enable_randomization_ )    return false;
//     if( do_continuous_randomization_!= collection.do_continuous_randomization_ )return false;
//     if( do_scattering_              != collection.do_scattering_ )           return false;
//     if( location_                   != collection.location_ )                return false;
//     if( density_correction_         != collection.density_correction_ )      return false;
//     if( do_exact_time_calculation_   != collection.do_exact_time_calculation_ )return false;
//     if( do_time_interpolation_      != collection.do_time_interpolation_ )   return false;
//     if( *time_particle_             != *collection.time_particle_ )          return false;
//
//     if( *decay_                     != *collection.decay_ )                  return false;
//
//     if( crosssections_.size()       != collection.crosssections_.size() )    return false;
//     if( bigLow_.size()              != collection.bigLow_.size() )           return false;
//     if( storeDif_.size()            != collection.storeDif_.size() )         return false;
//
//     for(unsigned int i =0; i<collection.bigLow_.size(); i++)
//     {
//         if( bigLow_.at(i) !=  collection.bigLow_.at(i) )        return false;
//     }
//
//     for(unsigned int i =0; i<collection.storeDif_.size(); i++)
//     {
//         if( storeDif_.at(i) !=  collection.storeDif_.at(i) )    return false;
//     }
//
//     for(unsigned int i =0; i<collection.crosssections_.size(); i++)
//     {
//         switch (collection.crosssections_.at(i)->GetType())
//         {
//             case ParticleType::Brems:
//                 if( *(Bremsstrahlung*)crosssections_.at(i) !=  *(Bremsstrahlung*)collection.crosssections_.at(i) )
//                 return false;
//                 break;
//             // case ParticleType::DeltaE:
//             //     if( *(Ionization*)crosssections_.at(i) != *(Ionization*)collection.crosssections_.at(i) ) return
//             false;
//             //     break;
//             // case ParticleType::EPair:
//             //     if( *(Epairproduction*)crosssections_.at(i) !=  *(Epairproduction*)collection.crosssections_.at(i)
//             ) return false;
//             //     break;
//             // case ParticleType::NuclInt:
//             //     if( *(Photonuclear*)crosssections_.at(i) !=  *(Photonuclear*)collection.crosssections_.at(i) )
//             return false;
//             //     break;
//             default:
//                 log_fatal("Unknown cross section");
//                 exit(1);
//         }
//     }
//
//     if( interpolant_ != NULL && collection.interpolant_ != NULL)
//     {
//         if( *interpolant_   != *collection.interpolant_)                                        return false;
//     }
//     else if( interpolant_ != collection.interpolant_)                                           return false;
//
//     if( interpolant_diff_ != NULL && collection.interpolant_diff_ != NULL)
//     {
//         if( *interpolant_diff_   != *collection.interpolant_diff_)                              return false;
//     }
//     else if( interpolant_diff_ != collection.interpolant_diff_)                                 return false;
//
//     if( interpol_prop_decay_ != NULL && collection.interpol_prop_decay_ != NULL)
//     {
//         if( *interpol_prop_decay_   != *collection.interpol_prop_decay_)                        return false;
//     }
//     else if( interpol_prop_decay_ != collection.interpol_prop_decay_)                           return false;
//
//     if( interpol_prop_decay_diff_ != NULL && collection.interpol_prop_decay_diff_ != NULL)
//     {
//         if( *interpol_prop_decay_diff_   != *collection.interpol_prop_decay_diff_)              return false;
//     }
//     else if( interpol_prop_decay_diff_ != collection.interpol_prop_decay_diff_)                 return false;
//
//     if( interpol_prop_interaction_ != NULL && collection.interpol_prop_interaction_ != NULL)
//     {
//         if( *interpol_prop_interaction_   != *collection.interpol_prop_interaction_)            return false;
//     }
//     else if( interpol_prop_interaction_ != collection.interpol_prop_interaction_)               return false;
//
//     if( interpol_prop_interaction_diff_ != NULL && collection.interpol_prop_interaction_diff_ != NULL)
//     {
//         if( *interpol_prop_interaction_diff_   != *collection.interpol_prop_interaction_diff_)  return false;
//     }
//     else if( interpol_prop_interaction_diff_ != collection.interpol_prop_interaction_diff_)     return false;
//
//     if( interpol_time_particle_diff_ != NULL && collection.interpol_time_particle_diff_ != NULL)
//     {
//         if( *interpol_time_particle_diff_   != *collection.interpol_time_particle_diff_)        return false;
//     }
//     else if( interpol_time_particle_diff_ != collection.interpol_time_particle_diff_)           return false;
//
//     if( interpol_time_particle_ != NULL && collection.interpol_time_particle_ != NULL)
//     {
//         if( *interpol_time_particle_   != *collection.interpol_time_particle_)                  return false;
//     }
//     else if( interpol_time_particle_ != collection.interpol_time_particle_)                     return false;
//
//     if( randomizer_ != NULL && collection.randomizer_ != NULL)
//     {
//         if( *randomizer_   != *collection.randomizer_)                                        return false;
//     }
//     else if( randomizer_ != collection.randomizer_)                                           return false;
//
//     if( scattering_ != NULL && collection.scattering_ != NULL)
//     {
//         if( *scattering_   != *collection.scattering_)                                        return false;
//     }
//     else if( scattering_ != collection.scattering_)                                           return false;
//
//     if( geometry_ != NULL && collection.geometry_ != NULL)
//     {
//         if( *geometry_   != *collection.geometry_)                                        return false;
//     }
//     else if( geometry_ != collection.geometry_)                                           return false;
//
//     //else
//     return true;
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// bool Collection::operator!=(const Collection &collection) const
// {
//     return !(*this == collection);
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// namespace PROPOSAL
// {
//
// ostream& operator<<(ostream& os, Collection const& collection)
// {
//
//     os<<"------------Collection( "<<&collection<<" )------------"<<endl;
//     os<<"------------------------------------------------------"<<endl;
//     os<<"Particle type:\t\t\t\t\t"<<collection.particle_->GetName()<<endl;
//     os<<"Order of interpolation:\t\t\t\t"<<collection.order_of_interpolation_<<endl;
//     os<<"Interpolation enabled:\t\t\t\t"<<collection.do_interpolation_<<endl;
//     os<<"Continuous randomization enabled:\t\t"<<collection.do_continuous_randomization_<<endl;
//     os<<"Moliere scattering enabled:\t\t\t"<<collection.do_scattering_<<endl;
//     os<<"Exact time calculation enabled:\t\t\t"<<collection.do_exact_time_calculation_<<endl;
//     os<<"Location (0=infront, 1=inside, 2=behind)\t"<<collection.location_<<endl;
//     os<<"Density correction factor:\t\t\t"<<collection.density_correction_<<endl;
//
//     os<<"\n"<<*collection.medium_<<"\n"<<endl;
//     os<<*collection.cut_settings_<<"\n"<<endl;
//
//     if(collection.geometry_ != NULL)
//         os<<*collection.geometry_<<"\n"<<endl;
//
//     os<<"------------------------------------------------------"<<endl;
//     os<<"------------------------------------------------------"<<endl;
//     return os;
// }
//
// }  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// void Collection::swap(Collection &collection)
// {
//     using std::swap;
//
//     PROPOSALParticle tmp_particle1(*collection.particle_);
//     PROPOSALParticle tmp_particle2(*particle_);
//
//     EnergyCutSettings tmp_cuts1(*collection.cut_settings_);
//     EnergyCutSettings tmp_cuts2(*cut_settings_);
//
//     vector<CrossSections*> tmp_cross1(collection.crosssections_);
//     vector<CrossSections*> tmp_cross2(crosssections_);
//
//     swap( order_of_interpolation_     , collection.order_of_interpolation_ );
//     swap( do_interpolation_           , collection.do_interpolation_ );
//     swap( lpm_effect_enabled_         , collection.lpm_effect_enabled_ );
//     swap( ini_                        , collection.ini_ );
//     swap( up_                         , collection.up_ );
//     swap( do_weighting_               , collection.do_weighting_ );
//     swap( weighting_order_            , collection.weighting_order_ );
//     swap( weighting_starts_at_        , collection.weighting_starts_at_ );
//     swap( do_continuous_randomization_, collection.do_continuous_randomization_ );
//     swap( enable_randomization_       , collection.enable_randomization_ );
//     swap( location_                   , collection.location_ );
//     swap( density_correction_         , collection.density_correction_ );
//     swap( do_time_interpolation_      , collection.do_time_interpolation_ );
//     swap( do_exact_time_calculation_   , collection.do_exact_time_calculation_ );
//
//
//     particle_->swap( *collection.particle_ );       //particle pointer swap
//     medium_->swap( *collection.medium_ );
//     integral_->swap( *collection.integral_ );
//     cut_settings_->swap( *collection.cut_settings_ );
//     crosssections_.swap(collection.crosssections_); //particle pointer swap
//     prop_decay_->swap( *collection.prop_decay_ );
//     prop_interaction_->swap( *collection.prop_interaction_ );
//     decay_->swap(*collection.decay_);               //particle pointer swap
//
//     time_particle_->swap(*collection.time_particle_ );
//
//     storeDif_.swap(collection.storeDif_);
//     bigLow_.swap(collection.bigLow_);
//
//     if( randomizer_ != NULL && collection.randomizer_ != NULL)
//     {
//         randomizer_->swap(*collection.randomizer_);
//     }
//     else if( randomizer_ == NULL && collection.randomizer_ != NULL)
//     {
//         randomizer_ = new ContinuousRandomization(*collection.randomizer_);
//         collection.randomizer_ = NULL;
//     }
//     else if( randomizer_ != NULL && collection.randomizer_ == NULL)
//     {
//         collection.randomizer_ = new ContinuousRandomization(*randomizer_);
//         randomizer_ = NULL;
//     }
//
//     if( scattering_ != NULL && collection.scattering_ != NULL)
//     {
//         scattering_->swap(*collection.scattering_);
//     }
//     else if( scattering_ == NULL && collection.scattering_ != NULL)
//     {
//         scattering_ = new Scattering(*collection.scattering_);
//         collection.scattering_ = NULL;
//     }
//     else if( scattering_ != NULL && collection.scattering_ == NULL)
//     {
//         collection.scattering_ = new Scattering(*scattering_);
//         scattering_ = NULL;
//     }
//
//     if( geometry_ != NULL && collection.geometry_ != NULL)
//     {
//         geometry_->swap(*collection.geometry_);
//     }
//     else if( geometry_ == NULL && collection.geometry_ != NULL)
//     {
//         // geometry_ = new Geometry(*collection.geometry_);
//         geometry_ = collection.geometry_->clone();
//         collection.geometry_ = NULL;
//     }
//     else if( geometry_ != NULL && collection.geometry_ == NULL)
//     {
//         // collection.geometry_ = new Geometry(*geometry_);
//         geometry_ = collection.geometry_->clone();
//         geometry_ = NULL;
//     }
//
//
//     // Set pointers again (to many swapping above....)
//     // SetParticle( new PROPOSALParticle(tmp_particle1) );
//     // collection.SetParticle( new PROPOSALParticle(tmp_particle2) );
//
//     SetCutSettings(  new EnergyCutSettings(tmp_cuts1) );
//     collection.SetCutSettings( new EnergyCutSettings(tmp_cuts2) );
//
//     SetCrosssections(  tmp_cross1 );
//     collection.SetCrosssections(  tmp_cross2 );
//
//     if( interpolant_ != NULL && collection.interpolant_ != NULL)
//     {
//         interpolant_->swap(*collection.interpolant_);
//     }
//     else if( interpolant_ == NULL && collection.interpolant_ != NULL)
//     {
//         interpolant_ = new Interpolant(*collection.interpolant_);
//         collection.interpolant_ = NULL;
//     }
//     else if( interpolant_ != NULL && collection.interpolant_ == NULL)
//     {
//         collection.interpolant_ = new Interpolant(*interpolant_);
//         interpolant_ = NULL;
//     }
//
//     if( interpolant_diff_ != NULL && collection.interpolant_diff_ != NULL)
//     {
//         interpolant_diff_->swap(*collection.interpolant_diff_);
//     }
//     else if( interpolant_diff_ == NULL && collection.interpolant_diff_ != NULL)
//     {
//         interpolant_diff_ = new Interpolant(*collection.interpolant_diff_);
//         collection.interpolant_diff_ = NULL;
//     }
//     else if( interpolant_diff_ != NULL && collection.interpolant_diff_ == NULL)
//     {
//         collection.interpolant_diff_ = new Interpolant(*interpolant_diff_);
//         interpolant_diff_ = NULL;
//     }
//
//     if( interpol_prop_decay_ != NULL && collection.interpol_prop_decay_ != NULL)
//     {
//         interpol_prop_decay_->swap(*collection.interpol_prop_decay_);
//     }
//     else if( interpol_prop_decay_ == NULL && collection.interpol_prop_decay_ != NULL)
//     {
//         interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_);
//         collection.interpol_prop_decay_ = NULL;
//     }
//     else if( interpol_prop_decay_ != NULL && collection.interpol_prop_decay_ == NULL)
//     {
//         collection.interpol_prop_decay_ = new Interpolant(*interpol_prop_decay_);
//         interpol_prop_decay_ = NULL;
//     }
//
//     if( interpol_prop_decay_diff_ != NULL && collection.interpol_prop_decay_diff_ != NULL)
//     {
//         interpol_prop_decay_diff_->swap(*collection.interpol_prop_decay_diff_);
//     }
//     else if( interpol_prop_decay_diff_ == NULL && collection.interpol_prop_decay_diff_ != NULL)
//     {
//         interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_);
//         collection.interpol_prop_decay_diff_ = NULL;
//     }
//     else if( interpol_prop_decay_diff_ != NULL && collection.interpol_prop_decay_diff_ == NULL)
//     {
//         collection.interpol_prop_decay_diff_ = new Interpolant(*interpol_prop_decay_diff_);
//         interpol_prop_decay_diff_ = NULL;
//     }
//
//     if( interpol_prop_interaction_ != NULL && collection.interpol_prop_interaction_ != NULL)
//     {
//         interpol_prop_interaction_->swap(*collection.interpol_prop_interaction_);
//     }
//     else if( interpol_prop_interaction_ == NULL && collection.interpol_prop_interaction_ != NULL)
//     {
//         interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_);
//         collection.interpol_prop_interaction_ = NULL;
//     }
//     else if( interpol_prop_interaction_ != NULL && collection.interpol_prop_interaction_ == NULL)
//     {
//         collection.interpol_prop_interaction_ = new Interpolant(*interpol_prop_interaction_);
//         interpol_prop_interaction_ = NULL;
//     }
//
//     if( interpol_prop_interaction_diff_ != NULL && collection.interpol_prop_interaction_diff_ != NULL)
//     {
//         interpol_prop_interaction_diff_->swap(*collection.interpol_prop_interaction_diff_);
//     }
//     else if( interpol_prop_interaction_diff_ == NULL && collection.interpol_prop_interaction_diff_ != NULL)
//     {
//         interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_);
//         collection.interpol_prop_interaction_diff_ = NULL;
//     }
//     else if( interpol_prop_interaction_diff_ != NULL && collection.interpol_prop_interaction_diff_ == NULL)
//     {
//         collection.interpol_prop_interaction_diff_ = new Interpolant(*interpol_prop_interaction_diff_);
//         interpol_prop_interaction_diff_ = NULL;
//     }
//
//     if( interpol_time_particle_ != NULL && collection.interpol_time_particle_ != NULL)
//     {
//         interpol_time_particle_->swap(*collection.interpol_time_particle_);
//     }
//     else if( interpol_time_particle_ == NULL && collection.interpol_time_particle_ != NULL)
//     {
//         interpol_time_particle_ = new Interpolant(*collection.interpol_time_particle_);
//         collection.interpol_time_particle_ = NULL;
//     }
//     else if( interpol_time_particle_ != NULL && collection.interpol_time_particle_ == NULL)
//     {
//         collection.interpol_time_particle_ = new Interpolant(*interpol_time_particle_);
//         interpol_time_particle_ = NULL;
//     }
//
//     if( interpol_time_particle_diff_ != NULL && collection.interpol_time_particle_diff_ != NULL)
//     {
//         interpol_time_particle_diff_->swap(*collection.interpol_time_particle_diff_);
//     }
//     else if( interpol_time_particle_diff_ == NULL && collection.interpol_time_particle_diff_ != NULL)
//     {
//         interpol_time_particle_diff_ = new Interpolant(*collection.interpol_time_particle_diff_);
//         collection.interpol_time_particle_diff_ = NULL;
//     }
//     else if( interpol_time_particle_diff_ != NULL && collection.interpol_time_particle_diff_ == NULL)
//     {
//         collection.interpol_time_particle_diff_ = new Interpolant(*interpol_time_particle_diff_);
//         interpol_time_particle_diff_ = NULL;
//     }
//
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Collection::FunctionToTimeIntegral(const PROPOSALParticle& particle, double energy)
{
    double aux;

    aux = FunctionToIntegral(particle, energy);
    aux *= particle.GetEnergy() / (particle.GetMomentum() * SPEED);
    return aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Collection::FunctionToPropIntegralDecay(const PROPOSALParticle& particle, double energy)
{
    double aux;
    double decay;

    aux = FunctionToIntegral(particle, energy);

    decay = MakeDecay(particle);

    log_debug(" + %f", particle.GetEnergy());

    return aux * decay;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Collection::FunctionToPropIntegralInteraction(const PROPOSALParticle& particle, double energy)
{
    double aux;
    double rate       = 0.0;
    double total_rate = 0.0;

    aux = FunctionToIntegral(particle, energy);

    for (unsigned int i = 0; i < crosssections_.size(); i++)
    {
        rate = crosssections_.at(i)->CalculatedNdx(particle);

        log_debug("Rate for %s = %f", crosssections_.at(i)->GetName().c_str(), rate);

        total_rate += rate;
    }
    return aux * total_rate;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Collection::FunctionToIntegral(const PROPOSALParticle& particle, double energy)
{
    double result;
    double aux;

    PROPOSALParticle temp_particle(particle);
    temp_particle.SetEnergy(energy);

    result = 0.0;

    for (unsigned int i = 0; i < crosssections_.size(); i++)
    {
        aux = crosssections_.at(i)->CalculatedEdx(temp_particle);
        result += aux;

        log_debug("energy %f , dE/dx = %f", particle.GetEnergy(), aux);
    }

    return -1 / result;
}
