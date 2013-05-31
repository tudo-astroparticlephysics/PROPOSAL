#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Constants.h"
#include <algorithm>
#include "boost/bind.hpp"
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double ContinuousRandomization::Randomize(double initial_energy, double final_energy, double rnd)
{
    return standard_normal_-> StandardNormalRandomNumber(
                rnd,
                final_energy,
                sqrt( DE2de(initial_energy, final_energy) ),
                particle_->GetLow(),
                initial_energy,
                false );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::EnableDE2dxInterpolation()
{
    if(do_dE2dx_Interpolation_)return;
    standard_normal_->EnableInterpolation();

    double energy = particle_->GetEnergy();
    dE2dx_interpolant_    =   new Interpolant(NUM2, particle_->GetLow(), BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2dxInterplant, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);

    do_dE2dx_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::EnableDE2deInterpolation()
{
    if(do_dE2de_Interpolation_)return;
    standard_normal_->EnableInterpolation();

    double energy = particle_->GetEnergy();

    dE2de_interpolant_       =   new Interpolant(NUM2, particle_->GetLow(), BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2deInterplant, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    dE2de_interpolant_diff_  =   new Interpolant(NUM2, particle_->GetLow(), BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2deInterplantDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);

    do_dE2de_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void ContinuousRandomization::DisableDE2dxInterpolation()
{
    delete dE2dx_interpolant_;
    do_dE2dx_Interpolation_=false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::DisableDE2deInterpolation()
{
    delete dE2de_interpolant_;
    delete dE2de_interpolant_diff_;
    do_dE2de_Interpolation_=false;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ContinuousRandomization::ContinuousRandomization()
    :cross_sections_            ( )
    ,do_dE2dx_Interpolation_    ( false )
    ,do_dE2de_Interpolation_    ( false )
    ,which_cross_               ( 0 )
    ,order_of_interpolation_    ( 5 )
{
    particle_           =   new Particle();
    medium_             =   new Medium();
    dE2dx_integral_     =   new Integral(IROMB, IMAXS, IPREC);
    dE2de_integral_     =   new Integral(IROMB, IMAXS, IPREC2);

    standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);

    dE2dx_interpolant_      =   NULL;
    dE2de_interpolant_      =   NULL;
    dE2de_interpolant_diff_ =   NULL;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ContinuousRandomization::ContinuousRandomization(Particle* particle, Medium* medium, std::vector<CrossSections*> cross_sections)
    :particle_                  ( particle )
    ,medium_                    ( medium )
    ,cross_sections_            ( cross_sections )
    ,do_dE2dx_Interpolation_    ( false )
    ,do_dE2de_Interpolation_    ( false )
    ,which_cross_               ( 0 )
    ,order_of_interpolation_    ( 5 )
{
    dE2dx_integral_     =   new Integral(IROMB, IMAXS, IPREC);
    dE2de_integral_     =   new Integral(IROMB, IMAXS, IPREC2);

    standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);

    dE2dx_interpolant_      =   NULL;
    dE2de_interpolant_      =   NULL;
    dE2de_interpolant_diff_ =   NULL;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ContinuousRandomization::ContinuousRandomization(const ContinuousRandomization &continuous_randomization)
    :particle_                  ( new Particle(*continuous_randomization.particle_) )
    ,medium_                    ( new Medium( *continuous_randomization.medium_) )
    ,do_dE2dx_Interpolation_    ( continuous_randomization.do_dE2dx_Interpolation_)
    ,do_dE2de_Interpolation_    ( continuous_randomization.do_dE2de_Interpolation_)
    ,standard_normal_           ( new StandardNormal(*continuous_randomization.standard_normal_) )
    ,dE2dx_integral_            ( new Integral(*continuous_randomization.dE2dx_integral_) )
    ,dE2de_integral_            ( new Integral(*continuous_randomization.dE2de_integral_) )
    ,which_cross_               ( continuous_randomization.which_cross_ )
    ,order_of_interpolation_    ( continuous_randomization.order_of_interpolation_ )
{
    cross_sections_.resize(continuous_randomization.cross_sections_.size());

    for(unsigned int i =0; i<continuous_randomization.cross_sections_.size(); i++)
    {
        if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            cross_sections_.at(i) = new Bremsstrahlung( *(Bremsstrahlung*)continuous_randomization.cross_sections_.at(i) );
        }
        else if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Ionization")==0)
        {
            cross_sections_.at(i) = new Ionization( *(Ionization*)continuous_randomization.cross_sections_.at(i) );
        }
        else if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            cross_sections_.at(i) = new Epairproduction( *(Epairproduction*)continuous_randomization.cross_sections_.at(i) );
        }
        else if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            cross_sections_.at(i) = new Photonuclear( *(Photonuclear*)continuous_randomization.cross_sections_.at(i) );
        }
        else
        {
            cout<<"In copy constructor of ContinuousRandomization: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    if(continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        dE2dx_interpolant_ = new Interpolant(*continuous_randomization.dE2dx_interpolant_) ;
    }
    else
    {
        dE2dx_interpolant_ = NULL;
    }

    if(continuous_randomization.dE2de_interpolant_ != NULL)
    {
        dE2de_interpolant_ = new Interpolant(*continuous_randomization.dE2de_interpolant_) ;
    }
    else
    {
        dE2de_interpolant_ = NULL;
    }

    if(continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        dE2de_interpolant_diff_ = new Interpolant(*continuous_randomization.dE2de_interpolant_diff_) ;
    }
    else
    {
        dE2de_interpolant_diff_ = NULL;
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ContinuousRandomization& ContinuousRandomization::operator=(const ContinuousRandomization &continuous_randomization)
{
    if (this != &continuous_randomization)
    {
      ContinuousRandomization tmp(continuous_randomization);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ContinuousRandomization::operator==(const ContinuousRandomization &continuous_randomization) const
{
    if( *particle_              != *continuous_randomization.particle_ )                return false;
    if( *medium_                != *continuous_randomization.medium_ )                  return false;
    if( *dE2dx_integral_        != *continuous_randomization.dE2dx_integral_ )          return false;
    if( *dE2de_integral_        != *continuous_randomization.dE2de_integral_ )          return false;
    if( do_dE2dx_Interpolation_ != continuous_randomization.do_dE2dx_Interpolation_ )   return false;
    if( do_dE2de_Interpolation_ != continuous_randomization.do_dE2de_Interpolation_ )   return false;
    if( *standard_normal_       != *continuous_randomization.standard_normal_ )         return false;
    if( which_cross_            != continuous_randomization.which_cross_)               return false;
    if( order_of_interpolation_ != continuous_randomization.order_of_interpolation_ )   return false;


    for(unsigned int i =0; i<continuous_randomization.cross_sections_.size(); i++)
    {
        if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            if( *(Bremsstrahlung*)cross_sections_.at(i) !=  *(Bremsstrahlung*)continuous_randomization.cross_sections_.at(i) ) return false;
        }
        else if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Ionization")==0)
        {
            if( *(Ionization*)cross_sections_.at(i) != *(Ionization*)continuous_randomization.cross_sections_.at(i) ) return false;
        }
        else if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            if( *(Epairproduction*)cross_sections_.at(i) !=  *(Epairproduction*)continuous_randomization.cross_sections_.at(i) ) return false;
        }
        else if(continuous_randomization.cross_sections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            if( *(Photonuclear*)cross_sections_.at(i) !=  *(Photonuclear*)continuous_randomization.cross_sections_.at(i) )  return false;
        }
        else
        {
            cout<<"In operator== of ContinuousRandomization: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    if( dE2dx_interpolant_ != NULL && continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        if( *dE2dx_interpolant_   != *continuous_randomization.dE2dx_interpolant_)                      return false;
    }
    else if( dE2dx_interpolant_ != continuous_randomization.dE2dx_interpolant_)                         return false;

    if( dE2de_interpolant_ != NULL && continuous_randomization.dE2de_interpolant_ != NULL)
    {
        if( *dE2de_interpolant_   != *continuous_randomization.dE2de_interpolant_)                      return false;
    }
    else if( dE2de_interpolant_ != continuous_randomization.dE2de_interpolant_)                         return false;

    if( dE2de_interpolant_diff_ != NULL && continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        if( *dE2de_interpolant_diff_   != *continuous_randomization.dE2de_interpolant_diff_)            return false;
    }
    else if( dE2de_interpolant_diff_ != continuous_randomization.dE2de_interpolant_diff_)               return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ContinuousRandomization::operator!=(const ContinuousRandomization &continuous_randomization) const
{
  return !(*this == continuous_randomization);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::swap(ContinuousRandomization &continuous_randomization)
{
    using std::swap;

    Particle tmp_particle1(*continuous_randomization.particle_);
    Particle tmp_particle2(*particle_);

    Medium tmp_medium1(*continuous_randomization.medium_);
    Medium tmp_medium2(*medium_);


    particle_->swap(*continuous_randomization.particle_);
    medium_->swap(*continuous_randomization.medium_);
    dE2dx_integral_->swap( *continuous_randomization.dE2dx_integral_ );
    dE2de_integral_->swap( *continuous_randomization.dE2de_integral_ );
    standard_normal_->swap( *continuous_randomization.standard_normal_ );

    cross_sections_.swap(continuous_randomization.cross_sections_);

    swap( which_cross_ , continuous_randomization.which_cross_ );
    swap( do_dE2dx_Interpolation_, continuous_randomization.do_dE2dx_Interpolation_);
    swap( do_dE2de_Interpolation_, continuous_randomization.do_dE2de_Interpolation_);
    swap( order_of_interpolation_, continuous_randomization.order_of_interpolation_);

    // Set pointers again (to many swapping above....)
    SetParticle( new Particle(tmp_particle1) );
    continuous_randomization.SetParticle( new Particle(tmp_particle2) );

    SetMedium( new Medium(tmp_medium1) );
    continuous_randomization.SetMedium( new Medium(tmp_medium2) );


    if( dE2dx_interpolant_ != NULL && continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        dE2dx_interpolant_->swap(*continuous_randomization.dE2dx_interpolant_);
    }
    else if( dE2dx_interpolant_ == NULL && continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        dE2dx_interpolant_ = new Interpolant(*continuous_randomization.dE2dx_interpolant_);
        continuous_randomization.dE2dx_interpolant_ = NULL;
    }
    else if( dE2dx_interpolant_ != NULL && continuous_randomization.dE2dx_interpolant_ == NULL)
    {
        continuous_randomization.dE2dx_interpolant_ = new Interpolant(*dE2dx_interpolant_);
        dE2dx_interpolant_ = NULL;
    }

    if( dE2de_interpolant_ != NULL && continuous_randomization.dE2de_interpolant_ != NULL)
    {
        dE2de_interpolant_->swap(*continuous_randomization.dE2de_interpolant_);
    }
    else if( dE2de_interpolant_ == NULL && continuous_randomization.dE2de_interpolant_ != NULL)
    {
        dE2de_interpolant_ = new Interpolant(*continuous_randomization.dE2de_interpolant_);
        continuous_randomization.dE2de_interpolant_ = NULL;
    }
    else if( dE2de_interpolant_ != NULL && continuous_randomization.dE2de_interpolant_ == NULL)
    {
        continuous_randomization.dE2de_interpolant_ = new Interpolant(*dE2de_interpolant_);
        dE2de_interpolant_ = NULL;
    }

    if( dE2de_interpolant_diff_ != NULL && continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        dE2de_interpolant_diff_->swap(*continuous_randomization.dE2de_interpolant_diff_);
    }
    else if( dE2de_interpolant_diff_ == NULL && continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        dE2de_interpolant_diff_ = new Interpolant(*continuous_randomization.dE2de_interpolant_diff_);
        continuous_randomization.dE2de_interpolant_diff_ = NULL;
    }
    else if( dE2de_interpolant_diff_ != NULL && continuous_randomization.dE2de_interpolant_diff_ == NULL)
    {
        continuous_randomization.dE2de_interpolant_diff_ = new Interpolant(*dE2de_interpolant_diff_);
        dE2de_interpolant_diff_ = NULL;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------private Memberfunctions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::DE2dx()
{
    if(do_dE2dx_Interpolation_)
    {
        return max( dE2dx_interpolant_->Interpolate(particle_->GetEnergy()) , 0.0 );
    }

    double sum = 0;
    double min = 0;

    for(unsigned int i = 0 ; i < cross_sections_.size() ; i++)
    {
        which_cross_ = i;

        for(int j=0 ; j < medium_->GetNumCompontents() ; j++)
        {
            cross_sections_.at(i)->SetIntegralLimits(j);

            if(cross_sections_.at(i)->GetName().compare("Bremsstrahlung")==0)
            {
                min =   0;
            }
            else
            {
                min =   cross_sections_.at(i)->GetVMin();
            }
            sum +=  dE2dx_integral_->Integrate (min, cross_sections_.at(i)->GetVUp(), boost::bind(&ContinuousRandomization::FunctionToDE2dxIntegral, this, _1) ,2);

            if(cross_sections_.at(i)->GetName().compare("Ionization")==0)
            {
                break;
            }
        }
    }

    return particle_->GetEnergy()*particle_->GetEnergy()*sum;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::DE2de( double ei, double ef )
{
    if( do_dE2de_Interpolation_ )
    {
        if( abs( ei-ef ) > abs(ei)*HALF_PRECISION )
        {
            double aux  =   dE2de_interpolant_->Interpolate( ei );
            double aux2 =   aux - dE2de_interpolant_->Interpolate( ef );

            if( abs(aux2) > abs(aux)*HALF_PRECISION )
            {
                return max(aux2, 0.0);
            }
        }

        else
        {
            return max( dE2de_interpolant_diff_->Interpolate( (ei+ef)/2 )*(ef-ei) , 0.0 );
        }

    }

    return dE2de_integral_->Integrate(ei, ef, boost::bind(&ContinuousRandomization::FunctionToDE2deIntegral, this, _1), 4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToBuildDE2dxInterplant(double energy)
{
    particle_->SetEnergy(energy);
    return DE2dx();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToBuildDE2deInterplant(double energy)
{
    return dE2de_integral_->Integrate(energy, particle_->GetLow(), boost::bind(&ContinuousRandomization::FunctionToDE2deIntegral, this, _1) , 4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToBuildDE2deInterplantDiff(double energy)
{
    return FunctionToDE2deIntegral(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToDE2dxIntegral(double v)
{
    return v*v*cross_sections_.at(which_cross_)->FunctionToDNdxIntegral(v);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToDE2deIntegral(double energy)
{

    double result;
    double aux;

    particle_->SetEnergy(energy);
    result  =    0;

    for(unsigned int i = 0 ; i<cross_sections_.size() ; i++)
    {
        aux     =   cross_sections_.at(i)->CalculatedEdx();
        result  +=  aux;
    }
    return -1/result*DE2dx();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------------Setter---------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::SetMedium(Medium* medium)
{
    medium_ = medium;
    for(unsigned int i = 0 ; i < cross_sections_.size() ; i++)
    {
        cross_sections_.at(i)->SetMedium(medium_);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::SetParticle(Particle* particle)
{
    particle_ = particle;
    for(unsigned int i = 0 ; i < cross_sections_.size() ; i++)
    {
        cross_sections_.at(i)->SetParticle(particle);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::SetCrosssections(
        std::vector<CrossSections*> crosssections) {
    cross_sections_ = crosssections;
}
