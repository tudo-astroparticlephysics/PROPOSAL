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
    :order_of_interpolation_    ( 5 )
    ,do_interpolation_          ( false )
    ,lpm_effect_enabled_        ( false )
    ,ini_                       ( 0 )
    ,debug_                     ( false )
    ,crosssections_             ( )
{

    interpolant_        = new Interpolant();
    interpolant_diff_   = new Interpolant();
    particle_           = new Particle();
    medium_             = new Medium();
    integral_           = new Integral();
    cut_settings_       = new EnergyCutSettings();

}
//----------------------------------------------------------------------------//

ProcessCollection::ProcessCollection(const ProcessCollection &collection)
    :order_of_interpolation_    ( collection.order_of_interpolation_ )
    ,do_interpolation_          ( collection.do_interpolation_ )
    ,lpm_effect_enabled_        ( collection.lpm_effect_enabled_ )
    ,ini_                       ( collection.ini_ )
    ,debug_                     ( collection.debug_ )
    ,interpolant_               ( new Interpolant (*collection.interpolant_) )
    ,interpolant_diff_          ( new Interpolant (*collection.interpolant_diff_) )
    ,particle_                  ( new Particle(*collection.particle_) )
    ,medium_                    ( new Medium( *collection.medium_) )
    ,integral_                  ( new Integral(*collection.integral_) )
    ,cut_settings_              ( new EnergyCutSettings(*collection.cut_settings_) )
{
    crosssections_.resize(collection.crosssections_.size());

    for(unsigned int i =0; i<collection.crosssections_.size(); i++)
    {
        if(collection.crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            crosssections_.at(i) = new Bremsstrahlung( *(Bremsstrahlung*)collection.crosssections_.at(i) );
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Ionization")==0)
        {
            crosssections_.at(i) = new Ionization( *(Ionization*)collection.crosssections_.at(i) );
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            crosssections_.at(i) = new Epairproduction( *(Epairproduction*)collection.crosssections_.at(i) );
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            crosssections_.at(i) = new Photonuclear( *(Photonuclear*)collection.crosssections_.at(i) );
        }
        else
        {
            cout<<"In copy constructor of ProcessCollection: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

}
//----------------------------------------------------------------------------//

ProcessCollection& ProcessCollection::operator=(const ProcessCollection &collection){
    if (this != &collection)
    {
      ProcessCollection tmp(collection);
      swap(tmp);
    }
    return *this;
}
//----------------------------------------------------------------------------//
bool ProcessCollection::operator==(const ProcessCollection &collection) const
{
    if( order_of_interpolation_    != collection.order_of_interpolation_ )  return false;
    if( do_interpolation_          != collection.do_interpolation_ )        return false;
    if( lpm_effect_enabled_        != collection.lpm_effect_enabled_ )      return false;
    if( ini_                       != collection.ini_ )                     return false;
    if( debug_                     != collection.debug_ )                   return false;
    if( *interpolant_              != *collection.interpolant_ )            return false;
    if( *interpolant_diff_         != *collection.interpolant_diff_ )       return false;
    if( *particle_                 != *collection.particle_ )               return false;
    if( *medium_                   != *collection.medium_ )                 return false;
    if( *integral_                 != *collection.integral_ )               return false;
    if( *cut_settings_             != *collection.cut_settings_ )           return false;

    if( crosssections_.size()      != collection.crosssections_.size() )    return false;

    for(unsigned int i =0; i<collection.crosssections_.size(); i++)
    {
        if(collection.crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            if( *(Bremsstrahlung*)crosssections_.at(i) !=  *(Bremsstrahlung*)collection.crosssections_.at(i) ) return false;
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Ionization")==0)
        {
            if( *(Ionization*)crosssections_.at(i) != *(Ionization*)collection.crosssections_.at(i) ) return false;

        }
        else if(collection.crosssections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            if( *(Epairproduction*)crosssections_.at(i) !=  *(Epairproduction*)collection.crosssections_.at(i) ) return false;
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            if( *(Photonuclear*)crosssections_.at(i) !=  *(Photonuclear*)collection.crosssections_.at(i) )  return false;
        }
        else
        {
            cout<<"In copy constructor of ProcessCollection: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    //else
    return true;
}
//----------------------------------------------------------------------------//
bool ProcessCollection::operator!=(const ProcessCollection &collection) const
{
    return !(*this == collection);
}
//----------------------------------------------------------------------------//
void ProcessCollection::swap(ProcessCollection &collection)
{
    using std::swap;

    swap( order_of_interpolation_    , collection.order_of_interpolation_ );
    swap( do_interpolation_          , collection.do_interpolation_ );
    swap( lpm_effect_enabled_        , collection.lpm_effect_enabled_ );
    swap( ini_                       , collection.ini_ );
    swap( debug_                     , collection.debug_ );

    interpolant_->swap( *collection.interpolant_ );
    interpolant_diff_->swap( *collection.interpolant_diff_ );
    particle_->swap( *collection.particle_ );
    medium_->swap( *collection.medium_ );
    integral_->swap( *collection.integral_ );
    cut_settings_->swap( *collection.cut_settings_ );
    crosssections_.swap(collection.crosssections_);

    for(unsigned int i = 0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->SetParticle(particle_);
        crosssections_.at(i)->SetMedium(medium_);
        crosssections_.at(i)->SetEnergyCutSettings(cut_settings_);

        collection.crosssections_.at(i)->SetParticle(collection.GetParticle());
        collection.crosssections_.at(i)->SetMedium(collection.GetMedium());
        collection.crosssections_.at(i)->SetEnergyCutSettings(collection.GetCutSettings());
    }

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
    crosssections_.at(0) = new Ionization(particle_, medium_, cut_settings_);
    crosssections_.at(1) = new Bremsstrahlung(particle_, medium_, cut_settings_);
    crosssections_.at(2) = new Photonuclear(particle_, medium_, cut_settings_);
    crosssections_.at(3) = new Epairproduction(particle_, medium_, cut_settings_);

    interpolant_        = new Interpolant();
    interpolant_diff_   = new Interpolant();

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


double ProcessCollection::GetDisplacement(double ei, double ef, double dist)
{
    if(do_interpolation_)
    {
        if(fabs(ei-ef) > fabs(ei)*HALF_PRECISION)
        {
            double aux;

            ini_    =   interpolant_->Interpolate(ei);
            aux     =   ini_ - interpolant_->Interpolate(ef);

            if(fabs(aux) > fabs(ini_)*HALF_PRECISION)
            {
                return max(aux, 0.0);
            }

        }

        ini_    =   0;
        return max((interpolant_diff_->Interpolate((ei + ef)/2))*(ef - ei), 0.0);

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
            aux     =   interpolant_->FindLimit(ini_-dist);

            if(fabs(aux) > fabs(ei)*HALF_PRECISION)
            {
                return min( max(aux, particle_->GetLow()), ei);
            }
        }

        return min( max(ei+dist/interpolant_diff_->Interpolate(ei + dist/(2*interpolant_diff_->Interpolate(ei))), particle_->GetLow()), ei);
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

void ProcessCollection::SetCrosssections(
		std::vector<CrossSections*> crosssections) {
	crosssections_ = crosssections;
}

void ProcessCollection::SetCutSettings(EnergyCutSettings* cutSettings) {
	cut_settings_ = cutSettings;
}

void ProcessCollection::SetDebug(bool debug) {
	debug_ = debug;
}

void ProcessCollection::SetDoInterpolation(bool doInterpolation) {
	do_interpolation_ = doInterpolation;
}

void ProcessCollection::SetIni(double ini) {
	ini_ = ini;
}

void ProcessCollection::SetIntegral(Integral* integral) {
	integral_ = integral;
}

void ProcessCollection::SetInterpolant(Interpolant* interpolant) {
	interpolant_ = interpolant;
}

void ProcessCollection::SetInterpolantDiff(
		Interpolant* interpolantDiff) {
	interpolant_diff_ = interpolantDiff;
}

void ProcessCollection::SetLpmEffectEnabled(bool lpmEffectEnabled) {
	lpm_effect_enabled_ = lpmEffectEnabled;
}

void ProcessCollection::SetMedium(Medium* medium) {
	medium_ = medium;
}

void ProcessCollection::SetOrderOfInterpolation(int orderOfInterpolation) {
	order_of_interpolation_ = orderOfInterpolation;
}

void ProcessCollection::SetParticle(Particle* particle) {
	particle_ = particle;
}

//----------------------------------------------------------------------------//

//Setter


//----------------------------------------------------------------------------//
//destructors

ProcessCollection::~ProcessCollection(){}
