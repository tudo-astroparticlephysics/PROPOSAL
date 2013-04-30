/*
 * ProcessCollection.cxx
 *
 *  Created on: 29.04.2013
 *      Author: koehne
 */

#include "PROPOSAL/ProcessCollection.h"
#include <cmath>
#include "boost/function.hpp"
#include "boost/bind.hpp"


using namespace std;


//constructors

ProcessCollection::ProcessCollection()
    :order_of_interpolation_(5)
    ,do_interpolation_(false)
    ,lpm_effect_enabled_(false)
    ,ini_(0)
    ,debug_(false)
{

}
//----------------------------------------------------------------------------//

ProcessCollection::ProcessCollection(const ProcessCollection &processcollection)
{
    *this = processcollection;
}
//----------------------------------------------------------------------------//

ProcessCollection& ProcessCollection::operator=(const ProcessCollection &processcollection){
    return *this;
}

//----------------------------------------------------------------------------//

ProcessCollection::ProcessCollection(Particle *particle, Medium *medium, EnergyCutSettings* cut_settings)
    :order_of_interpolation_(5)
    ,do_interpolation_(false)
    ,lpm_effect_enabled_(false)
    ,ini_(0)
    ,debug_(false)
{

    particle_       =   particle;
    medium_         =   medium;
    cut_settings_   =   cut_settings;

    integral_       =   new Integral(IROMB, IMAXS, IPREC2);

    crosssections_.resize(4);
    crosssections_.at(1) = new Ionization(particle_, medium_, cut_settings_);
    crosssections_.at(2) = new Bremsstrahlung(particle_, medium_, cut_settings_);
    crosssections_.at(3) = new Photonuclear(particle_, medium_, cut_settings_);
    crosssections_.at(4) = new Epairproduction(particle_, medium_, cut_settings_);

}

//Memberfunctions

//----------------------------------------------------------------------------//

double ProcessCollection::FunctionToBuildInterpolant(double energy)
{
    return integral_->Integrate(energy, particle_->GetLow(), boost::bind(&ProcessCollection::FunctionToIntegral, this, _1),4);
}

//----------------------------------------------------------------------------//

double ProcessCollection::FunctionToBuildInterpolantDiff(double energy)
{
    return FunctionToIntegral(energy);
}


//----------------------------------------------------------------------------//

double ProcessCollection::FunctionToIntegral(double energy)
{

    double result;
    double aux;

    particle_->SetEnergy(energy);
    result  =    0;

    if(debug_)
    {
        cout<<" * "<<particle_->GetEnergy();
    }

    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        aux     =   crosssections_.at(i)->CalculatedEdx();
        result  +=  aux;

        if(debug_)
        {
            cout<<" \t "<<aux;
        }
    }

    return -1/result;
}

//----------------------------------------------------------------------------//


double ProcessCollection::GetDx(double ei, double ef, double dist)
{
    if(do_interpolation_)
    {
        if(fabs(ei-ef) > fabs(ei)*HALF_PRECISION)
        {
            double aux;

            ini_    =   interpolant_->interpolate(ei);
            aux     =   ini_ - interpolant_->interpolate(ef);

            if(fabs(aux) > fabs(ini_)*HALF_PRECISION)
            {
                return max(aux, 0.0);
            }

        }

        ini_    =   0;
        return max((interpolant_diff_->interpolate((ei + ef)/2))*(ef - ei), 0.0);

    }
    else
    {
        return integral_->IntegrateWithLog(ei, ef, boost::bind(&ProcessCollection::FunctionToIntegral, this, _1), -dist);
    }

}

//----------------------------------------------------------------------------//

double ProcessCollection::GetEf(double ei, double dist)
{

    if(do_interpolation_)
    {
        if(ini_ != 0)
        {
            double aux;
            aux     =   interpolant_->findLimit(ini_-dist);

            if(fabs(aux) > fabs(ei)*HALF_PRECISION)
            {
                return min( max(aux, particle_->GetLow()), ei);
            }
        }

        return min( max(ei+dist/interpolant_diff_->interpolate(ei + dist/(2*interpolant_diff_->interpolate(ei))), particle_->GetLow()), ei);
    }
    else
    {
        return integral_->GetUpperLimit();
    }

}


//----------------------------------------------------------------------------//


void ProcessCollection::EnableInterpolation()
{
    if(do_interpolation_)return;

    double energy = particle_->GetEnergy();

    interpolant_        =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::FunctionToBuildInterpolant, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpolant_diff_   =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::FunctionToBuildInterpolantDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);

    do_interpolation_=true;
}

//----------------------------------------------------------------------------//

void ProcessCollection::DisableInterpolation()
{
    do_interpolation_  =   false;

}

//----------------------------------------------------------------------------//

void ProcessCollection::EnableLpmEffect()
{
    lpm_effect_enabled_ =true;
    for(unsigned int i=0;i<crosssections_.size();i++){
        crosssections_.at(i)->EnableLpmEffect(lpm_effect_enabled_);
    }
}

//----------------------------------------------------------------------------//

void ProcessCollection::DisableLpmEffect()
{
    lpm_effect_enabled_ =false;
    for(unsigned int i=0;i<crosssections_.size();i++){
        crosssections_.at(i)->EnableLpmEffect(lpm_effect_enabled_);
    }
}

//----------------------------------------------------------------------------//

//Setter


//----------------------------------------------------------------------------//
//destructors

ProcessCollection::~ProcessCollection(){}
