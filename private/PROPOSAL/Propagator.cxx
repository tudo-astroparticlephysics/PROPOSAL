/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
#include <utility>
#include <boost/program_options.hpp>

using namespace std;

namespace po	= boost::program_options;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


std::pair<double,double> Propagator::CalculateEnergyTillStochastic( double initial_energy )
{
    double rndd    =-  log(RandomDouble());
    double rndi    =-  log(RandomDouble());

    double rndiMin = 0;
    double rnddMin = 0;

    pair<double,double> final;

    //solving the tracking integral
    if(particle_->GetLifetime() < 0)
    {
        rnddMin =   0;
    }
    else
    {
        rnddMin =   collection_->CalculateTrackingIntegal(initial_energy, rndd, false)/density_correction_;
    }

    rndiMin =   collection_->CalculateTrackingIntegal(initial_energy, rndi, true);

    //evaluating the energy loss
    if(rndd >= rnddMin || rnddMin<=0)
    {
        final.second =   particle_->GetLow();
    }
    else
    {
        final.second =   collection_->CalculateFinalEnergy( initial_energy, rndd*density_correction_, false );
    }

    if(rndi >= rndiMin || rndiMin <= 0)
    {
        final.first =   particle_->GetLow();
    }
    else
    {
        final.first =   collection_->CalculateFinalEnergy( initial_energy, rndi, true );
    }

    return final;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::Propagate( double distance )
{
    bool    flag;
    double  displacement;

    double  initial_energy  =   particle_->GetEnergy();
    double  final_energy    =   particle_->GetEnergy();

    pair<double,string> decay;
    pair<double,string> energy_loss;

    //first: final energy befor first interaction second: energy at which the
    // particle decay
    //first and second are compared to decide if interaction happens or decay
    pair<double,double> energy_till_stochastic_;


    if(distance < 0)
    {
        distance   =   0;
    }

    if(initial_energy <= particle_->GetLow() || distance==0)
    {
        flag    =   false;
    }
    else
    {
        flag    =   true;
    }

    while(flag)
    {
        energy_till_stochastic_ = CalculateEnergyTillStochastic( initial_energy );

        if(energy_till_stochastic_.first > energy_till_stochastic_.second)
        {
            particle_interaction_   =   true;
            final_energy            =   energy_till_stochastic_.first;
        }
        else
        {
            particle_interaction_   =   false;
            final_energy            =   energy_till_stochastic_.second;

        }
        cout<<"efi "<<energy_till_stochastic_.first<<"\t"<<energy_till_stochastic_.second<<"\t";

        cout<<final_energy<<"\t";

        //Calculate the displacement according to initial energy initial_energy and final_energy
        displacement  =   collection_->CalculateDisplacement(
                    initial_energy,
                    final_energy,
                    density_correction_*(distance - particle_->GetPropagatedDistance())) / density_correction_;
        cout<<particle_->GetT()<<"\t";
        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if( displacement > distance - particle_->GetPropagatedDistance() )
        {
            displacement  =   distance - particle_->GetPropagatedDistance();

            final_energy  =   collection_->CalculateFinalEnergy(initial_energy, density_correction_*displacement);

        }

        //Advance the Particle according to the displacement
        //Initial energy and final energy are needed if Molier Scattering is enabled
        AdvanceParticle(displacement, initial_energy, final_energy);

        if(abs(distance - particle_->GetPropagatedDistance()) < abs(distance)*COMPUTER_PRECISION)
        {
            particle_->SetPropagatedDistance( distance );  // computer precision control
        }

//        if(contiCorr)
//        {
//            if(ef!= particle_->low)
//            {
//                ef  =   StandardN->sndrn(RandomDouble(), ef, sqrt(E2Loss_->e2le->dE2de(ei, ef)), particle_->low, ei, false);
//            }

//        }

        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if( final_energy == particle_->GetLow() || particle_->GetPropagatedDistance() == distance)
        {
            break;
        }

        //Set the particle energy to the current energy before making
        //stochatic losses or decay
        particle_->SetEnergy( final_energy );

        if(particle_interaction_)
        {
            energy_loss     =   collection_->MakeStochasticLoss();
            final_energy    -=  energy_loss.first;

            cout<<energy_loss.first<<"\t"<<energy_loss.second<<endl;
        }
        else
        {
            decay           =   collection_->MakeDecay();
            final_energy    =   0;

            cout<<decay.first<<"\t"<<decay.second<<endl;
        }

        //break if the lower limit of particle energy is reached
        if(final_energy <= particle_->GetLow())
        {

            break;
        }

        //Next round: update the inital energy
        initial_energy  =   final_energy;

    }

//    if(sdec)
//    {
//        if(particle_->r!=r && ef!=0 && particle_->l>=0)
//        {
//            particle_->setEnergy(particle_->m);

//            particle_->t    +=  -particle_->l*log(RandomDouble());

//            if(particle_->type==2)
//            {
//                aux =   cros->get_decay()->e(RandomDouble(), 0.5, RandomDouble(), o);
//            }
//            else
//            {
//                aux =   cros->get_decay()->e(RandomDouble(), 0.5, 0.5, o);
//            }

//            ef  =   0;

//            o->output(1, cros->get_decay()->get_out(), aux, ef);
//        }
//    }

    particle_->SetEnergy(final_energy);

    //Particle reached the border, final energy is returned
    if(particle_->GetPropagatedDistance()==distance)
    {
        return final_energy;
    }
    //The particle stopped/decayed, the propageted distance is return with a minus sign
    else
    {
        return -particle_->GetPropagatedDistance();
    }
    //Should never be here
    return 0;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::CalculateParticleTime(double ei, double ef)
{
    if(do_time_interpolation_)
    {
        if(abs(ei-ef) > abs(ei)*HALF_PRECISION)
        {
            double aux  =   interpol_time_particle_->Interpolate(ei);
            double aux2 =   aux - interpol_time_particle_->Interpolate(ef);

            if(abs(aux2) > abs(aux)*HALF_PRECISION)
            {
                return aux2;
            }
        }

        return interpol_time_particle_diff_->Interpolate( (ei+ef)/2 )*(ef-ei);
    }
    else
    {
        return time_particle_->Integrate(ei, ef, boost::bind(&Propagator::FunctionToTimeIntegral, this, _1),4);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::AdvanceParticle(double dr, double ei, double ef)
{

    double dist = particle_->GetPropagatedDistance();
    double time = particle_->GetT();
    double x    = particle_->GetX();
    double y    = particle_->GetY();
    double z    = particle_->GetZ();

    dist   +=  dr;

    if(do_exact_time_calulation_)
    {
        time   +=  CalculateParticleTime(ei, ef)/density_correction_;
    }
    else
    {
        time   +=  dr/SPEED;
    }


//    if(propagate_->get_molieScat())
//    Implement the Molie Scattering here see PROPOSALParticle::advance of old version

//    else
//    {
    x   +=  particle_->GetSinTheta() * particle_->GetCosPhi() * dr;
    y   +=  particle_->GetSinTheta() * particle_->GetSinPhi() * dr;
    z   +=  particle_->GetCosTheta() * dr;

    particle_->SetPropagatedDistance(dist);
    particle_->SetT(time);
    particle_->SetX(x);
    particle_->SetY(y);
    particle_->SetZ(z);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::EnableParticleTimeInterpolation()
{
    if(do_time_interpolation_)return;

    double energy = particle_->GetEnergy();

    interpol_time_particle_         =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&Propagator::InterpolTimeParticle, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpol_time_particle_diff_    =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&Propagator::InterpolTimeParticleDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);
    do_time_interpolation_ =true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::DisableParticleTimeInterpolation()
{
    delete interpol_time_particle_;
    delete interpol_time_particle_diff_;

    interpol_time_particle_         = NULL;
    interpol_time_particle_diff_    = NULL;

    do_time_interpolation_ =false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::EnableInterpolation()
{
    collection_->EnableInterpolation();
    EnableParticleTimeInterpolation();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::DisableInterpolation()
{
    collection_->DisableInterpolation();
    DisableParticleTimeInterpolation();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Set and validate options--------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


boost::program_options::options_description Propagator::CreateOptions()
{

    po::options_description general("General options");
    general.add_options()
        ("help,h",		"shows this message")
        ("version,v",	"shows the version of the program");


    po::options_description propagator("Propagator options");
    propagator.add_options()
        ("propagator.interpol_time",        po::value<bool>(&do_time_interpolation_)->implicit_value(false),    "Enables interpolation for particle time calculation")
        ("propagator.exact_time",           po::value<bool>(&do_exact_time_calulation_)->implicit_value(false), "Do exact particle time calculations")
        ("propagator.density_correction",   po::value<double>(&density_correction_)->default_value(1),          "Density correction factor")
        ("propagator.interpol_order",       po::value<int>(&order_of_interpolation_)->default_value(5),         "number of interpolation points");


    po::options_description all("All options");
        all.add(general);
        all.add(propagator);

    for(unsigned int i =0 ; i < collection_->GetCrosssections().size(); i++)
    {
        all.add(collection_->GetCrosssections().at(i)->CreateOptions());
    }

    return all;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::Setup(int argc, char** argv)
{
    po::options_description all = CreateOptions();

    //parse cmd line
    po::variables_map vm;
    po::store( po::command_line_parser(argc, argv).options(all).run(), vm);

    //print help message if wanted
    if(vm.count("help")) {
        std::cout<< all;
        exit(1);
    }
    //notifies globalVar
    try {
        //set the variables
        vm.notify();
    }
    catch (po::invalid_command_line_syntax &e) {
        std::cerr<<"Error: "<<e.what()<<"\n";
        exit(1);
    }

    for(unsigned int i = 0 ; i < collection_->GetCrosssections().size() ; i++)
    {
        collection_->GetCrosssections().at(i)->ValidateOptions();
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator()
    :order_of_interpolation_    ( 5 )
    ,debug_                     ( false )
    ,particle_interaction_      ( false )
    ,density_correction_        ( 1. )
    ,do_time_interpolation_     ( false )
    ,do_exact_time_calulation_  ( false )
{
    particle_              = new Particle("mu",0,0,0,0,0,1e6,0);
    time_particle_         = new Integral();

    interpol_time_particle_         = NULL;
    interpol_time_particle_diff_    = NULL;

    InitDefaultCollection();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator(const Propagator &propagator)
    :order_of_interpolation_    ( propagator.order_of_interpolation_ )
    ,debug_                     ( propagator.debug_ )
    ,particle_interaction_      ( propagator.particle_interaction_ )
    ,density_correction_        ( propagator.density_correction_ )
    ,do_time_interpolation_     ( propagator.do_time_interpolation_ )
    ,do_exact_time_calulation_  ( propagator.do_exact_time_calulation_ )
    ,particle_                  ( propagator.particle_ )
    ,collection_                ( new ProcessCollection(*propagator.collection_) )
    ,time_particle_             ( new Integral(*propagator.time_particle_) )

{
    if(propagator.interpol_time_particle_ != NULL)
    {
        interpol_time_particle_ = new Interpolant(*propagator.interpol_time_particle_) ;
    }
    else
    {
        interpol_time_particle_ = NULL;
    }

    if(propagator.interpol_time_particle_diff_ != NULL)
    {
        interpol_time_particle_diff_ = new Interpolant(*propagator.interpol_time_particle_diff_) ;
    }
    else
    {
        interpol_time_particle_diff_ = NULL;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator& Propagator::operator=(const Propagator &propagator)
{
    if (this != &propagator)
    {
      Propagator tmp(propagator);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Propagator::operator==(const Propagator &propagator) const
{
    if( order_of_interpolation_   != propagator.order_of_interpolation_ ) return false;
    if( debug_                    != propagator.debug_ )                  return false;
    if( particle_                 != propagator.particle_ )               return false;
    if( particle_interaction_     != propagator.particle_interaction_ )   return false;
    if( density_correction_       != propagator.density_correction_ )     return false;
    if( do_time_interpolation_    != propagator.do_time_interpolation_ )  return false;
    if( do_exact_time_calulation_ != propagator.do_exact_time_calulation_ )return false;

    if( *collection_              != *propagator.collection_ )            return false;
    if( *time_particle_           != *propagator.time_particle_ )         return false;

    if( interpol_time_particle_diff_ != NULL && propagator.interpol_time_particle_diff_ != NULL)
    {
        if( *interpol_time_particle_diff_   != *propagator.interpol_time_particle_diff_)        return false;
    }
    else if( interpol_time_particle_diff_ != propagator.interpol_time_particle_diff_)           return false;

    if( interpol_time_particle_ != NULL && propagator.interpol_time_particle_ != NULL)
    {
        if( *interpol_time_particle_   != *propagator.interpol_time_particle_)                  return false;
    }
    else if( interpol_time_particle_ != propagator.interpol_time_particle_)                     return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Propagator::operator!=(const Propagator &propagator) const
{
    return !(*this == propagator);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::swap(Propagator &propagator)
{
    using std::swap;

    swap( order_of_interpolation_   ,   propagator.order_of_interpolation_ );
    swap( debug_                    ,   propagator.debug_);
    swap( particle_interaction_     ,   propagator.particle_interaction_);
    swap( density_correction_       ,   propagator.density_correction_ );
    swap( do_time_interpolation_    ,   propagator.do_time_interpolation_ );
    swap( do_exact_time_calulation_ ,   propagator.do_exact_time_calulation_ );


    particle_->swap( *propagator.particle_ );
    collection_->swap( *propagator.collection_ );
    time_particle_->swap(*propagator.time_particle_ );

    if( interpol_time_particle_ != NULL && propagator.interpol_time_particle_ != NULL)
    {
        interpol_time_particle_->swap(*propagator.interpol_time_particle_);
    }
    else if( interpol_time_particle_ == NULL && propagator.interpol_time_particle_ != NULL)
    {
        interpol_time_particle_ = new Interpolant(*propagator.interpol_time_particle_);
        propagator.interpol_time_particle_ = NULL;
    }
    else if( interpol_time_particle_ != NULL && propagator.interpol_time_particle_ == NULL)
    {
        propagator.interpol_time_particle_ = new Interpolant(*interpol_time_particle_);
        interpol_time_particle_ = NULL;
    }

    if( interpol_time_particle_diff_ != NULL && propagator.interpol_time_particle_diff_ != NULL)
    {
        interpol_time_particle_diff_->swap(*propagator.interpol_time_particle_diff_);
    }
    else if( interpol_time_particle_diff_ == NULL && propagator.interpol_time_particle_diff_ != NULL)
    {
        interpol_time_particle_diff_ = new Interpolant(*propagator.interpol_time_particle_diff_);
        propagator.interpol_time_particle_diff_ = NULL;
    }
    else if( interpol_time_particle_diff_ != NULL && propagator.interpol_time_particle_diff_ == NULL)
    {
        propagator.interpol_time_particle_diff_ = new Interpolant(*interpol_time_particle_diff_);
        interpol_time_particle_diff_ = NULL;
    }


}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::InitDefaultCollection()
{
    Medium* med             = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(500,-1);
    collection_             = new ProcessCollection(particle_ , med, cuts);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::InterpolTimeParticle(double energy)
{
    return FunctionToTimeIntegral(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::InterpolTimeParticleDiff(double energy)
{
    return time_particle_->Integrate(energy, particle_->GetLow(), boost::bind(&Propagator::FunctionToTimeIntegral, this, _1),4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::FunctionToTimeIntegral(double energy)
{
    double aux;

    aux     =   collection_->FunctionToIntegral(energy);
    aux     *=  particle_->GetEnergy()/(particle_->GetMomentum()*SPEED);
    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Propagator::~Propagator(){}





