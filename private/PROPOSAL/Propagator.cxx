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


using namespace std;


Propagator::Propagator()
    :debug_                 ( false )
    ,particle_interaction_  ( false )
    ,rho_                   ( 1. )
{
    particle_              = new Particle("mu",0,0,0,0,0,1e6,0);
    InitDefaultCollection();
}
//----------------------------------------------------------------------------//

Propagator::Propagator(const Propagator &propagator)
    :debug_                 ( propagator.debug_ )
    ,particle_interaction_  ( propagator.particle_interaction_ )
    ,rho_                   ( propagator.rho_ )
    ,particle_              ( propagator.particle_ )
    ,collection_            ( new ProcessCollection(*propagator.collection_) )
{

}

//----------------------------------------------------------------------------//

Propagator& Propagator::operator=(const Propagator &propagator){
    if (this != &propagator)
    {
      Propagator tmp(propagator);
      swap(tmp);
    }
    return *this;
}
//----------------------------------------------------------------------------//
bool Propagator::operator==(const Propagator &propagator) const
{
    if( debug_                  != propagator.debug_ )                  return false;
    if( particle_               != propagator.particle_ )               return false;
    if( particle_interaction_   != propagator.particle_interaction_ )   return false;
    if( rho_                    != propagator.rho_ )                    return false;
    if( *collection_            != *propagator.collection_ )            return false;
    //else
    return true;
}
//----------------------------------------------------------------------------//
bool Propagator::operator!=(const Propagator &propagator) const
{
    return !(*this == propagator);
}


//----------------------------------------------------------------------------//
void Propagator::swap(Propagator &propagator)
{
    using std::swap;

    swap( debug_                 ,   propagator.debug_);
    swap( particle_interaction_  ,   propagator.particle_interaction_);
    swap( rho_                   ,   propagator.rho_ );
    particle_->swap( *propagator.particle_ );
    collection_->swap( *propagator.collection_ );

}
//----------------------------------------------------------------------------//

void Propagator::InitDefaultCollection()
{
    Medium* med             = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(500,-1);
    collection_             = new ProcessCollection(particle_ , med, cuts);

}

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
        rnddMin =   collection_->CalculateTrackingIntegal(initial_energy, rndd, false)/rho_;
    }

    rndiMin =   collection_->CalculateTrackingIntegal(initial_energy, rndi, true);

    //evaluating the energy loss
    if(rndd >= rnddMin || rnddMin<=0)
    {
        final.second =   particle_->GetLow();
    }
    else
    {
        final.second =   collection_->CalculateFinalEnergy( initial_energy, rndd*rho_, false );
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
double Propagator::Propagate( double distance )
{
    bool    flag;
    double  displacement;

    double  initial_energy  =   particle_->GetEnergy();
    double  final_energy    =   particle_->GetEnergy();

    pair<double,string> decay;
    pair<double,string> energy_loss;

    //first: final energy befor first interaction second: decay
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
        displacement  =   collection_->CalculateDisplacement(initial_energy, final_energy, rho_*(distance - particle_->GetPropagatedDistance())) / rho_;

        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if( displacement > distance - particle_->GetPropagatedDistance() )
        {
            displacement  =   distance - particle_->GetPropagatedDistance();

            final_energy  =   collection_->CalculateFinalEnergy(initial_energy, rho_*displacement);

        }

        //Advance the Particle according to the displacement
        //Initial energy and final energy are needed if Molier Scattering is enabled
        collection_->AdvanceParticle(displacement, initial_energy, final_energy);

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
Propagator::~Propagator(){}





