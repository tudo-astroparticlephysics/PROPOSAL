#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Constants.h"
#include <algorithm>
#include "boost/bind.hpp"

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ContinuousRandomization::ContinuousRandomization()
{
    particle_           =   new Particle();
    medium_             =   new Medium();
    standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ContinuousRandomization::ContinuousRandomization(const ContinuousRandomization &continuous_randomization)
{
    particle_   = new Particle( *continuous_randomization.particle_ );
    medium_     = new Medium( *continuous_randomization.medium_ );
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
    if( *particle_  != *continuous_randomization.particle_ )             return false;
    if( *medium_    != *continuous_randomization.medium_ )               return false;

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

    particle_->swap(*continuous_randomization.particle_);
    medium_->swap(*continuous_randomization.medium_);
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
    return cross_sections_.at(which_cross_)->FunctionToDNdxIntegral(v);
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


void ContinuousRandomization::SetParticle(Particle *particle){
    particle_ = particle;
}

void ContinuousRandomization::SetMedium(Medium *medium){
    medium_ = medium;
}

void ContinuousRandomization::SetCrosssections(
        std::vector<CrossSections*> crosssections) {
    cross_sections_ = crosssections;
}
