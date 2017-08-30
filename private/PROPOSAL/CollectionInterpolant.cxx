
#include <boost/bind.hpp>

#include "PROPOSAL/CollectionInterpolant.h"
#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
// #include "PROPOSAL/Constants.h"

using namespace PROPOSAL;
using namespace std;

CollectionInterpolant::CollectionInterpolant()
    : Collection()
    , initialized_interpolation_(false)
    , up_(false)
    , big_low_interatction_(0)
    , big_low_decay_(0)
    , store_dif_interaction_(0)
    , store_dif_decay_(0)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpol_time_particle_(NULL)
    , interpol_time_particle_diff_(NULL)
    , interpol_prop_decay_(NULL)
    , interpol_prop_decay_diff_(NULL)
    , interpol_prop_interaction_(NULL)
    , interpol_prop_interaction_diff_(NULL)
{}

CollectionInterpolant::CollectionInterpolant(const Medium& medium,
                                             const Geometry& geometry,
                                             const EnergyCutSettings& cut_settings,
                                             const CollectionDef& def)
    :Collection(medium, geometry, cut_settings, def)
    , initialized_interpolation_(false)
    , up_(false)
    , big_low_interatction_(0)
    , big_low_decay_(0)
    , store_dif_interaction_(0)
    , store_dif_decay_(0)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpol_time_particle_(NULL)
    , interpol_time_particle_diff_(NULL)
    , interpol_prop_decay_(NULL)
    , interpol_prop_decay_diff_(NULL)
    , interpol_prop_interaction_(NULL)
    , interpol_prop_interaction_diff_(NULL)
{
}

CollectionInterpolant::~CollectionInterpolant()
{
    delete interpolant_;
    delete interpolant_diff_;
    delete interpol_time_particle_;
    delete interpol_time_particle_diff_;
    delete interpol_prop_decay_;
    delete interpol_prop_decay_diff_;
    delete interpol_prop_interaction_;
    delete interpol_prop_interaction_diff_;
}

CollectionInterpolant::CollectionInterpolant(const CollectionInterpolant& collection)
    :Collection(collection)
     ,initialized_interpolation_(collection.initialized_interpolation_)
     ,up_(collection.up_)
     ,big_low_interatction_(collection.big_low_interatction_)
     ,big_low_decay_(collection.big_low_decay_)
{
    // Only copy interpolants if there were already initialized in collections
    if (collection.initialized_interpolation_)
    {
        interpolant_ = new Interpolant(*collection.interpolant_);
        interpolant_diff_ = new Interpolant(*collection.interpolant_diff_);
        interpol_time_particle_ = new Interpolant(*collection.interpol_time_particle_);
        interpol_time_particle_diff_ = new Interpolant(*collection.interpol_time_particle_diff_);
        interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_);
        interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_);
        interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_);
        interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_);
    }

    // if (collection.interpolant_ != NULL)
    //     interpolant_ = new Interpolant(*collection.interpolant_);
    // if (collection.interpolant_diff_ != NULL)
    //     interpolant_diff_ = new Interpolant(*collection.interpolant_diff_);
    // if (collection.interpol_time_particle_ != NULL)
    //     interpol_time_particle_ = new Interpolant(*collection.interpol_time_particle_);
    // if (collection.interpol_time_particle_diff_ != NULL)
    //     interpol_time_particle_diff_ = new Interpolant(*collection.interpol_time_particle_diff_);
    // if (collection.interpol_prop_decay_ != NULL)
    //     interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_);
    // if (collection.interpol_prop_decay_diff_ != NULL)
    //     interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_);
    // if (collection.interpol_prop_interaction_ != NULL)
    //     interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_);
    // if (collection.interpol_prop_interaction_diff_ != NULL)
    //     interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_);
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::Propagate(PROPOSALParticle& particle, double distance)
{
    // Initialize Interpolants once first
    if (initialized_interpolation_ == false)
    {
        InitInterpolation(particle, collection_def_.path_to_tables, collection_def_.raw);
    }

    return Collection::Propagate(particle, distance);
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::CalculateDisplacement(const PROPOSALParticle& particle, double ei, double ef, double dist)
{
    (void)particle;
    (void)dist;

    if (fabs(ei - ef) > fabs(ei) * HALF_PRECISION)
    {
        double aux;

        ini_ = interpolant_->Interpolate(ei);
        aux  = ini_ - interpolant_->Interpolate(ef);

        if (fabs(aux) > fabs(ini_) * HALF_PRECISION)
        {
            return std::max(aux, 0.0);
        }
    }

    ini_ = 0;

    return std::max((interpolant_diff_->Interpolate((ei + ef) / 2)) * (ef - ei), 0.0);
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::CalculateFinalEnergy(const PROPOSALParticle& particle, double ei, double dist)
{
    if (ini_ != 0)
    {
        double aux;
        aux = interpolant_->FindLimit(ini_ - dist);

        if (fabs(aux) > fabs(ei) * HALF_PRECISION)
        {
            return std::min(std::max(aux, particle.GetLow()), ei);
        }
    }

    return std::min(
        std::max(ei + dist / interpolant_diff_->Interpolate(ei + dist / (2 * interpolant_diff_->Interpolate(ei))),
                 particle.GetLow()),
        ei);
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::CalculateFinalEnergy(const PROPOSALParticle& particle,
                                                double ei,
                                                double rnd,
                                                bool particle_interaction)
{
    if( particle_interaction )
    {
        if( abs(rnd) > abs(store_dif_interaction_)*HALF_PRECISION)
        {
            double aux;

            if(up_&&particle_interaction)
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_-rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_-rnd);
                }
            }
            else
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_+rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_+rnd);
                }
            }
            if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
            {
                return std::min(std::max(aux, particle.GetLow()), ei);
            }
        }
    }
    else
    {
        if(abs(rnd) > abs(store_dif_decay_)*HALF_PRECISION)
        {
            double aux;

            if(up_&&particle_interaction)
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_-rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_-rnd);
                }
            }
            else
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_+rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_+rnd);
                }
            }

            if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
            {
                return std::min(std::max(aux, particle.GetLow()), ei);
            }
        }
    }

    if(particle_interaction)
    {
        return std::min(std::max(ei + rnd/interpol_prop_interaction_diff_->Interpolate(ei + rnd/(2*interpol_prop_interaction_diff_->Interpolate(ei))), particle.GetLow()), ei);
    }
    else
    {
        return std::min(std::max(ei + rnd/interpol_prop_decay_diff_->Interpolate(ei + rnd/(2*interpol_prop_decay_diff_->Interpolate(ei))), particle.GetLow()), ei);
    }
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::CalculateTrackingIntegal(const PROPOSALParticle& particle,
                                                    double initial_energy,
                                                    double rnd,
                                                    bool particle_interaction)
{
    (void) particle;
    (void) rnd;

    if(particle_interaction)
    {
        store_dif_interaction_ =   interpol_prop_interaction_->Interpolate(initial_energy);
    }
    else
    {
        store_dif_decay_ =   interpol_prop_decay_->Interpolate(initial_energy);
    }

    if(up_&&particle_interaction)
    {
        if(particle_interaction)
        {
            return std::max(store_dif_interaction_, 0.0);
        }
        else
        {
            return std::max(store_dif_decay_, 0.0);
        }
    }
    else
    {
        if(particle_interaction)
        {
            return std::max(big_low_interatction_-store_dif_interaction_, 0.0);
        }
        else
        {
            return std::max(big_low_decay_-store_dif_decay_, 0.0);
        }
    }
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::CalculateParticleTime(const PROPOSALParticle& particle, double ei, double ef)
{
    (void) particle;

    if (abs(ei - ef) > abs(ei) * HALF_PRECISION)
    {
        double aux  = interpol_time_particle_->Interpolate(ei);
        double aux2 = aux - interpol_time_particle_->Interpolate(ef);

        if (abs(aux2) > abs(aux) * HALF_PRECISION)
        {
            return aux2;
        }
    }

    return interpol_time_particle_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
}

// ------------------------------------------------------------------------- //
void CollectionInterpolant::InitInterpolation(const PROPOSALParticle& particle, std::string filepath, bool raw)
{
    Integral prop_interaction_(IROMB, IMAXS, IPREC2);

    for(unsigned int i =0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->EnableDEdxInterpolation(particle, filepath,raw);
    }

    for(unsigned int i =0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->EnableDNdxInterpolation(particle, filepath,raw);
    }

    bool reading_worked =   true;
    bool storing_failed =   false;

    double a = abs(-prop_interaction_.Integrate(particle.GetLow(), particle.GetLow()*10, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),4));
    double b = abs(-prop_interaction_.Integrate(BIGENERGY, BIGENERGY/10, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),4));

    if( a < b)
    {
        up_  =   true;
    }
    else
    {
        up_  =   false;
    }

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables

    if(!filepath.empty())
    {
        std::stringstream filename;
        filename<<filepath<<"/Collection"
                <<"_"<<particle.GetName()
                <<"_mass_"<<particle.GetMass()
                <<"_charge_"<<particle.GetCharge()
                <<"_lifetime_"<<particle.GetLifetime()
                <<"_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_"<<cut_settings_.GetEcut()
                <<"_"<<cut_settings_.GetVcut();

        for(unsigned int i =0; i<crosssections_.size(); i++)
        {
            switch (crosssections_.at(i)->GetType())
            {
                case ParticleType::Brems:
                    filename << "_b"
                        << "_" << crosssections_.at(i)->GetParametrization()
                        << "_" << crosssections_.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::DeltaE:
                    filename << "_i";
                    break;
                case ParticleType::EPair:
                    filename << "_e"
                        << "_" << crosssections_.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::NuclInt:
                    filename << "_p"
                        << "_" << crosssections_.at(i)->GetParametrization();
                    break;
                default:
                    log_fatal("Unknown cross section");
                    exit(1);
            }
            filename<< "_" << crosssections_.at(i)->GetMultiplier()
                    << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetEcut()
                    << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetVcut();
        }

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_info("ProcessCollection parametrisation tables will be read from file:\t%s",filename.str().c_str());
            std::ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), std::ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            interpolant_        =   new Interpolant();
            interpolant_diff_   =   new Interpolant();

            reading_worked = interpolant_->Load(input,raw);
            reading_worked = interpolant_diff_->Load(input,raw);


            interpol_prop_decay_            =   new Interpolant();
            interpol_prop_decay_diff_       =   new Interpolant();
            interpol_prop_interaction_      =   new Interpolant();
            interpol_prop_interaction_diff_ =   new Interpolant();

            reading_worked = interpol_prop_decay_->Load(input,raw);
            reading_worked = interpol_prop_decay_diff_->Load(input,raw);
            reading_worked = interpol_prop_interaction_->Load(input,raw);
            reading_worked = interpol_prop_interaction_diff_->Load(input,raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {

            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("ProcessCollection parametrisation tables will be saved to file:\t%s",filename.str().c_str());

            if(abs(-prop_interaction_.Integrate(particle.GetLow(), particle.GetLow()*10, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),4))
                    < abs(-prop_interaction_.Integrate(BIGENERGY, BIGENERGY/10, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),4)))
            {
                up_  =   true;
            }
            else
            {
                up_  =   false;
            }

            std::ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), std::ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }

            if(output.good())
            {
                output.precision(16);

                double order_of_interpolation = collection_def_.order_of_interpolation;

                interpolant_        =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToBuildInterpolant, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
                interpolant_diff_   =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToIntegral, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);

                interpol_prop_decay_            =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::InterpolPropDecay, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
                interpol_prop_decay_diff_       =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralDecay, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
                interpol_prop_interaction_      =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::InterpolPropInteraction, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
                interpol_prop_interaction_diff_ =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);

                interpolant_->Save(output,raw);
                interpolant_diff_->Save(output,raw);
                interpol_prop_decay_->Save(output,raw);
                interpol_prop_decay_diff_->Save(output,raw);
                interpol_prop_interaction_->Save(output,raw);
                interpol_prop_interaction_diff_->Save(output,raw);

            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(filepath.empty() || storing_failed)
    {
        log_info("ProcessCollection parametrisation tables will be stored in memory!");

        double order_of_interpolation = collection_def_.order_of_interpolation;

        interpolant_        =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToBuildInterpolant, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
        interpolant_diff_   =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToIntegral, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);


        interpol_prop_decay_            =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::InterpolPropDecay, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
        interpol_prop_decay_diff_       =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralDecay, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
        interpol_prop_interaction_      =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::InterpolPropInteraction, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
        interpol_prop_interaction_diff_ =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);

    }

    if(collection_def_.do_continuous_randomization_)
    {
        randomizer_->EnableDE2dxInterpolation(particle, crosssections_, filepath ,raw);
        randomizer_->EnableDE2deInterpolation(particle, crosssections_, filepath,raw);
    }

    if(collection_def_.do_exact_time_calculation)
    {
        InitTimeInterpolation(particle, filepath,raw);
    }

    if(collection_def_.do_scattering)
    {
        scattering_->EnableInterpolation(particle, crosssections_, filepath);
    }

    big_low_decay_        = interpol_prop_decay_->Interpolate(particle.GetLow());
    big_low_interatction_ = interpol_prop_interaction_->Interpolate(particle.GetLow());

    initialized_interpolation_ = true;
}

// ------------------------------------------------------------------------- //
void CollectionInterpolant::InitTimeInterpolation(const PROPOSALParticle& particle, std::string filepath, bool raw)
{
    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables
    std::string particle_name = particle.GetName();

    if(!filepath.empty())
    {
        std::stringstream filename;
        filename << filepath << "/Time"
                 << "_" << particle_name << "_mass_" << particle.GetMass() << "_charge_" << particle.GetCharge()
                 << "_lifetime_" << particle.GetLifetime() << "_" << medium_->GetName() << "_"
                 << medium_->GetMassDensity() << "_" << cut_settings_.GetEcut() << "_" << cut_settings_.GetVcut();

        for(unsigned int i =0; i<crosssections_.size(); i++)
        {
            switch (crosssections_.at(i)->GetType())
            {
                case ParticleType::Brems:
                    filename << "_b"
                        << "_" << crosssections_.at(i)->GetParametrization()
                        << "_" << crosssections_.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::DeltaE:
                    filename << "_i";
                    break;
                case ParticleType::EPair:
                    filename << "_e"
                        << "_" << crosssections_.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::NuclInt:
                    filename << "_p"
                        << "_" << crosssections_.at(i)->GetParametrization();
                    break;
                default:
                    log_fatal("Unknown cross section");
                    exit(1);
            }
            filename<< "_" << crosssections_.at(i)->GetMultiplier()
                    << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetEcut()
                    << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetVcut();

        }

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Particle time parametrisation tables will be read from file:\t%s",filename.str().c_str());
            std::ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), std::ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            interpol_time_particle_         = new Interpolant();
            interpol_time_particle_diff_    = new Interpolant();

            reading_worked = interpol_time_particle_->Load(input,raw);
            reading_worked = interpol_time_particle_diff_->Load(input,raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Particle time parametrisation tables will be saved to file:\t%s",filename.str().c_str());

            std::ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), std::ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }

            if(output.good())
            {
                output.precision(16);

                double order_of_interpolation = collection_def_.order_of_interpolation;

                interpol_time_particle_         =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToTimeIntegral, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
                interpol_time_particle_diff_    =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);

                interpol_time_particle_->Save(output,raw);
                interpol_time_particle_diff_->Save(output,raw);

            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(filepath.empty() || storing_failed)
    {
        double order_of_interpolation = collection_def_.order_of_interpolation;

        interpol_time_particle_         =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToTimeIntegral, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);
        interpol_time_particle_diff_    =   new Interpolant(NUM3, particle.GetLow(), BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1), order_of_interpolation, false, false, true, order_of_interpolation, false, false, false);

    }
}

// ------------------------------------------------------------------------- //
// Helper functions to init interpolation
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CollectionInterpolant::FunctionToBuildInterpolant(const PROPOSALParticle& particle, double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    return integral.Integrate(energy, particle.GetLow(), boost::bind(&CollectionInterpolant::FunctionToIntegral, this, boost::cref(particle), _1),4);
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::InterpolPropDecay(const PROPOSALParticle& particle, double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    return -integral.Integrate(energy, BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralDecay, this, boost::cref(particle), _1),4);
}


// ------------------------------------------------------------------------- //
double CollectionInterpolant::InterpolPropInteraction(const PROPOSALParticle& particle, double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    if(up_)
    {
        return integral.Integrate(energy, particle.GetLow(), boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),4);
    }
    else
    {
        return -integral.Integrate(energy, BIGENERGY, boost::bind(&CollectionInterpolant::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),4);
    }
}

// ------------------------------------------------------------------------- //
double CollectionInterpolant::InterpolTimeParticleDiff(const PROPOSALParticle& particle, double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    return integral.Integrate(energy, particle.GetLow(), boost::bind(&CollectionInterpolant::FunctionToTimeIntegral, this, boost::cref(particle), _1),4);
}
