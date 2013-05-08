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
    ,rates_                 ( )
    ,total_rate_            ( 0 )
    ,particle_              ( )
    ,collection_            ( )
{

}
//----------------------------------------------------------------------------//

Propagator::Propagator(const Propagator &propagator)
    :debug_                 ( propagator.debug_ )
    ,particle_interaction_  ( propagator.particle_interaction_ )
    ,rho_                   ( propagator.rho_ )
    ,do_weighting_          ( propagator.do_weighting_ )
    ,weighting_order_       ( propagator.weighting_order_ )
    ,weighting_starts_at_   ( propagator.weighting_starts_at_ )
    ,rates_                 ( propagator.rates_ )
    ,total_rate_            ( propagator.total_rate_ )
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
    if( do_weighting_           != propagator.do_weighting_ )           return false;
    if( weighting_order_        != propagator.weighting_order_ )        return false;
    if( weighting_starts_at_    != propagator.weighting_starts_at_ )    return false;
    if( rates_.size()           != propagator.rates_.size() )           return false;
    if( total_rate_             != propagator.total_rate_ )             return false;

    for( unsigned int i = 0; i < propagator.rates_.size(); i++)
    {
        if( rates_.at(i)        != propagator.rates_.at(i) )            return false;
    }


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
    swap( do_weighting_          ,   propagator.do_weighting_ );
    swap( weighting_order_       ,   propagator.weighting_order_ );
    swap( weighting_starts_at_   ,   propagator.weighting_starts_at_ );
    swap( total_rate_            ,   propagator.total_rate_ );

    rates_.swap(propagator.rates_);
    particle_->swap( *propagator.particle_ );
    collection_->swap( *propagator.collection_ );

}
//----------------------------------------------------------------------------//

double Propagator::Propagate(double distance, double energy)
{
    //int wint;
    bool  flag;
    double ei, ef=0, efd, displacement, efi;//, aux=0;
    double rndd, rndi,rnddMin , rnd2, rnd3,  rndiMin, rnd1, rndTot;

    ei  =   energy;
    ef  =   ei;

    if(distance < 0)
    {
        distance   =   0;
    }


    if(energy <= particle_->GetEnergy() || distance==0)
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
            efd =   collection_->CalculateFinalEnergy(ei, rndd*rho, false);
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

            ef  =   collection_->GetEf(ei, rho_*displacement);

            if(debug_)
            {
                cerr<<"lost "<<ei - ef<<" MeV  ef = "<<ef<<" MeV"<<endl;
            }

            if(debug_)
            {
                cerr<<"5. calculating the local time ...  "<<endl;
            }
        }

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

        rnd2    =   RandomDouble();
        rnd3    =   RandomDouble();

        particle_->SetEnergy( ef );

        if(particle_interaction_)
        {
            rnd1    =   RandomDouble();

            if(do_weighting_)
            {
                if(particle_->GetPropagationDistance() > weighting_starts_at_)
                {
                    double exp      =   abs(weighting_order_);
                    double power    =   pow(rnd2, exp);

                    if(weighting_order_>0)
                    {
                        rnd2    =   1 - power*rnd2;
                    }
                    else
                    {
                        rnd2    =   power*rnd2;
                    }

                    weighting_order_        =   (1 + exp)*power;
                    weighting_starts_at_    =   particle_->GetPropagationDistance();
                    do_weighting_           =   false;
                }
            }

            rates_.resize(collection_->GetCrosssections().size());

            for(unsigned int i = 0 ; i < collection_->GetCrosssections().size(); i++)
            {
                rates_.at(i) =      collection_->GetCrosssections().at(i)->CalculatedNdx( rnd2 );
                total_rate_ +=  rates_.at(i);
            }

            rndTot = total_rate_*rnd1;

            for(unsigned int i = 0 ; i < rates_.size(); i++)
            {
                if(rates_.at(i) > rndTot)
                {
                    aux     =   collection_->GetCrosssections().at(i)->e(rnd3);
                    ef      -=  aux;
                    wint    =   2;
                }
            }

//            if(debug_)
//            {
//                decayS  =   cros->get_decay()->decay();
//            }


//            if(debug_)
//            {
//                cerr<<" . rnd1 = "<<o->f(rnd1)<<" rnd2 = "<<o->f(rnd2)<<
//                    " rnd3 = "<<o->f(rnd3)<<" decay = "<<o->f(decayS)<<endl;
//            }

//            if(debug_)
//            {
//                cerr<<" . ioniz = "<<o->f(ionizS)<<" brems = "<<o->f(bremsS)<<
//                    " photo = "<<o->f(photoS)<<" epair = "<<o->f(epairS);
//            }

//            for(unsigned int i = 0 ; i < rates_.size(); i++)

//            if(rates_.at(i)>rndTot)
//            {

//                aux     =   collection_->GetCrosssections().at(i)->
//                ef      -=  aux;
//                wint    =   2;
//            }
//            else if(ionizS+bremsS>rndTot)
//            {
//                aux     =   cros->get_bremsstrahlung()->get_Stochastic()->e(rnd3);
//                ef      -=  aux;
//                wint    =   3;
//            }
//            else if(ionizS+bremsS+photoS>rndTot)
//            {
//                aux     =   cros->get_photonuclear()->get_Stochastic()->e(rnd3);
//                ef      -=  aux;
//                wint    =   4;
//            }
//            else if(ionizS+bremsS+photoS+epairS>rndTot)
//            {
//                aux     =   cros->get_epairproduction()->stochastic_->e(rnd3);
//                ef      -=  aux;
//                wint    =   5;
//            }
//            else  // due to the parameterization of the cross section cutoffs
//            {
//                ei  =   ef;
//                continue;
//            }

//        }

//        else
//        {
//            if(particle_->type==2)
//            {
//                aux     =   cros->get_decay()->e(rnd2, rnd3, RandomDouble(), o);
//                ef      =   0;
//                wint    =   1;
//            }
//            else
//            {
//                aux     =   cros->get_decay()->e(rnd2, rnd3, 0.5, o);
//                ef      =   0;
//                wint    =   1;
//            }
//        }

//        if(wint==1)
//        {
//            o->output(wint, cros->get_decay()->get_out(), aux, ef);
//        }
//        else
//        {
//            o->output(wint, medium_->get_E(cros->get_component()), aux, ef);
//        }



//        if(ef<=particle_->low)
//        {

//            break;
//        }

//        ei  =   ef;

//    }

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

//    particle_->setEnergy(ef);  // to remember const state of the particle
//    o->HIST =   -1;            // to make sure user resets particle properties

//    if(particle_->r==r)
//    {
//        if(debug_)
//        {
//            cerr<<"PROPOSALParticle reached the border with energy ef = "<<o->f(ef)<<" MeV";
//        }

//        return ef;
//    }
//    else
//    {
//        if(debug_)
//        {
//            if(particle_->l<0)
//            {
//                cerr<<"PROPOSALParticle stopped at rf = "<<o->f(particle_->r)<<" cm";
//            }
//            else
//            {
//                cerr<<"PROPOSALParticle disappeared at rf = "<<o->f(particle_->r)<<" cm";
//            }
//        }

        return -particle_->GetPropagationDistance();
    }
    return 0;
}

//----------------------------------------------------------------------------//
Propagator::~Propagator(){}





