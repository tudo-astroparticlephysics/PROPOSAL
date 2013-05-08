/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Constants.h"
#include <cmath>


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
double Propagator::Propagate(double distance, double energy)
{
    //int wint;
    bool  flag;
    double ei, ef=0, efd, displacement, efi;//, aux=0;
    double energy_loss;
    double rndd, rndi,rnddMin , rndiMin;

    ei  =   energy;
    ef  =   ei;

    if(distance < 0)
    {
        distance   =   0;
    }


    if(energy <= particle_->GetLow() || distance==0)
    {
        flag    =   false;
    }
    else
    {
        flag    =   true;
    }

    if(debug_)
    {

        cerr<<"\nPropagating "<<particle_->GetName()<<" of energy "<<ei<<" MeV to a distance of "<<distance<<" cm"<<endl;
    }
    while(flag)
    {

        rndd    =-  log(RandomDouble());
        rndi    =-  log(RandomDouble());
        if(debug_)
        {
            cerr<<"1. solving the tracking integral ...  rndd = "<<rndd<<"  rndi = "<<rndi<<" ...  "<<endl;
        }

        if(particle_->GetLifetime() < 0)
        {
            rnddMin =   0;
        }
        else
        {
            rnddMin =   collection_->CalculateTrackingIntegal(ei, rndd, false)/rho_;
        }

        if(debug_)
        {
            cerr<<" \t \t \t rnddMin = "<<rnddMin<<" (d)  "<<endl;
        }

        rndiMin =   collection_->CalculateTrackingIntegal(ei, rndi, true);

        if(debug_)
        {
            cerr<<"rndiMin = "<<rndiMin<<" (i)"<<endl;
        }

        if(debug_)
        {
            cerr<<"2. evaluating the energy loss ...  "<<endl;
        }

        if(rndd >= rnddMin || rnddMin<=0)
        {
            efd =   particle_->GetLow();
        }
        else
        {
            efd =   collection_->CalculateFinalEnergy(ei, rndd*rho_, false);
        }

        if(debug_)
        {
            cerr<<"efd = "<<efd<<" MeV  "<<endl;
        }

        if(rndi >= rndiMin || rndiMin <= 0)
        {
            efi =   particle_->GetLow();
        }
        else
        {
            efi =   collection_->CalculateFinalEnergy(ei, rndi, true);
        }

        if(debug_)
        {
            cerr<<"efi = "<<efi<<" MeV ...  "<<endl;
        }
        particle_interaction_    =   (efi > efd);
        cout<<efi<<"\t"<<efd<<"\t";
        if(particle_interaction_)
        {
            ef  =   efi;
        }
        else
        {
            ef  =   efd;
        }

        if(debug_)
        {
            cerr<<" \t \t \t lost "<<ei-ef<<" MeV  ef = "<<ef<<" MeV"<<endl;
        }

        if(debug_)
        {
            cerr<<"3. calculating the displacement ...  "<<endl;
        }

        displacement  =   collection_->CalculateDisplacement(ei, ef, rho_*(distance - particle_->GetPropagationDistance())) / rho_;

        if(debug_)
        {
            cerr<<"displacement = "<<displacement<<" cm"<<endl;
        }

        if( displacement < distance - particle_->GetPropagationDistance() )
        {
            if(debug_)
            {
                cerr<<"4. calculating the local time ...  "<<endl;
            }
        }
        else
        {
            displacement  =   distance - particle_->GetPropagationDistance();

            if(debug_)
            {
                cerr<<"4. getting the const energy ...  "<<endl;
            }

            ef  =   collection_->CalculateFinalEnergy(ei, rho_*displacement);
            if(debug_)
            {
                cerr<<"lost "<<ei - ef<<" MeV  ef = "<<ef<<" MeV"<<endl;
            }

            if(debug_)
            {
                cerr<<"5. calculating the local time ...  "<<endl;
            }
        }
        cout<<ef<<endl;

//        if(recc)
//        {
//            o->output(0, "a"+particle_->name, ei, displacement);
//        }

//        particle_->advance(displacement, ei, ef);

        if(debug_)
        {
            cerr<<"t = "<<particle_->GetT()<<" s"<<endl;
        }

        if(abs(distance - particle_->GetPropagationDistance()) < abs(distance)*COMPUTER_PRECISION)
        {
            particle_->SetPropagationDistance( distance );  // computer precision control
        }

//        if(contiCorr)
//        {
//            if(ef!= particle_->low)
//            {
//                ef  =   StandardN->sndrn(RandomDouble(), ef, sqrt(E2Loss_->e2le->dE2de(ei, ef)), particle_->low, ei, false);
//            }

//        }

//        if(recc)
//        {
//            o->output(0, "conti", ei-ef, -displacement);
//        }

        if( ef == particle_->GetLow() || particle_->GetPropagationDistance() == distance)
        {
            break;
        }

        if(debug_)
        {
            cerr<<"5. choosing the cross section ..."<<endl;
        }

        energy_loss = collection_->MakeStochasticLoss(particle_interaction_ , ef);
        if(energy_loss<0)cout<<"UUUU"<<endl;
        if(ef <= particle_->GetLow())
        {

            break;
        }

        ei  =   ef;

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

    particle_->SetEnergy(ef);  // to remember const state of the particle
//    o->HIST =   -1;            // to make sure user resets particle properties

    if(particle_->GetPropagationDistance()==distance)
    {
        if(debug_)
        {
            cerr<<"PROPOSALParticle reached the border with energy ef = "<<ef<<" MeV";
        }

        return ef;
    }
    else
    {
        if(debug_)
        {
            if(particle_->GetLifetime()<0)
            {
                cerr<<"PROPOSALParticle stopped at rf = "<<particle_->GetPropagationDistance()<<" cm";
            }
            else
            {
                cerr<<"PROPOSALParticle disappeared at rf = "<<particle_->GetPropagationDistance()<<" cm";
            }
        }

        return -particle_->GetPropagationDistance();
    }
    return 0;
}

//----------------------------------------------------------------------------//
Propagator::~Propagator(){}





