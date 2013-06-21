/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Scattering.h"

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Scattering::Scattering( )
    :x0_(0)
    ,do_interpolation_(false)
    ,order_of_interpolation_(5)
    ,integral_(  new Integral(IROMB, IMAXS, IPREC2) )
    ,interpolant_(NULL)
    ,interpolant_diff_(NULL)
    ,particle_( new Particle("mu") )
{
    crosssections_.push_back(new Ionization());
    crosssections_.push_back(new Bremsstrahlung());
    crosssections_.push_back(new Photonuclear());
    crosssections_.push_back(new Epairproduction());

    SetParticle(particle_);

    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        if( crosssections_.at(i)->GetName().compare("Bremsstrahlung") ==0){
            Bremsstrahlung* tmp = (Bremsstrahlung*)crosssections_.at(i);
            x0_ = tmp->CalculateScatteringX0();
        }
    }
}

//double Scattering::cutoff   =   1;


Scattering::Scattering(std::vector<CrossSections*> crosssections)
    :x0_(0)
    ,do_interpolation_(false)
    ,order_of_interpolation_(5)
    ,integral_( new Integral(IROMB, IMAXS, IPREC2) )
    ,interpolant_(NULL)
    ,interpolant_diff_(NULL)
    ,particle_(NULL)
    ,crosssections_( crosssections)
{
    SetParticle(crosssections.at(0)->GetParticle());
    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        if( crosssections_.at(i)->GetName().compare("Bremsstrahlung") ==0){
            Bremsstrahlung* tmp = (Bremsstrahlung*)crosssections_.at(i);
            x0_ = tmp->CalculateScatteringX0();
        }
    }
}

Scattering::Scattering(const Scattering &scattering)
    :x0_(scattering.x0_)
    ,do_interpolation_(scattering.do_interpolation_)
    ,particle_(scattering.particle_)
    ,crosssections_(scattering.crosssections_)
{
    if(scattering.interpolant_ != NULL)
    {
        interpolant_ = new Interpolant(*scattering.interpolant_) ;
    }
    else
    {
        interpolant_ = NULL;
    }

    if(scattering.interpolant_diff_ != NULL)
    {
        interpolant_diff_ = new Interpolant(*scattering.interpolant_diff_) ;
    }
    else
    {
        interpolant_diff_ = NULL;
    }

    integral_ = new Integral(* scattering.integral_);
};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Scattering& Scattering::operator=(const Scattering &scattering){
    if (this != &scattering)
    {
      Scattering tmp(scattering);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator==(const Scattering &scattering) const
{
    if(*particle_ != *(scattering.particle_))return false;

    for(unsigned int i =2; i<scattering.crosssections_.size(); i++)
    {
        if(scattering.crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            if( *(Bremsstrahlung*)crosssections_.at(i) !=  *(Bremsstrahlung*)scattering.crosssections_.at(i) ) return false;
        }
        else if(scattering.crosssections_.at(i)->GetName().compare("Ionization")==0)
        {
            if( *(Ionization*)crosssections_.at(i) != *(Ionization*)scattering.crosssections_.at(i) ) return false;
        }
        else if(scattering.crosssections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            if( *(Epairproduction*)crosssections_.at(i) !=  *(Epairproduction*)scattering.crosssections_.at(i) ) return false;
        }
        else if(scattering.crosssections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            if( *(Photonuclear*)crosssections_.at(i) !=  *(Photonuclear*)scattering.crosssections_.at(i) )  return false;
        }
        else
        {
            cout<<"In copy constructor of Scattering: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    if( interpolant_ != NULL && scattering.interpolant_ != NULL)
    {
        if( *interpolant_   != *scattering.interpolant_)                                        return false;
    }

    if( interpolant_diff_ != NULL && scattering.interpolant_diff_ != NULL)
    {
        if( *interpolant_diff_   != *scattering.interpolant_diff_)                                        return false;
    }

    if( integral_ != NULL && scattering.integral_ != NULL)
    {
        if( *integral_   != *scattering.integral_)                                        return false;
    }

    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator!=(const Scattering &scattering) const {
  return !(*this == scattering);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Scattering::swap(Scattering &scattering)
{
    using std::swap;

    Particle tmp_particle1(*scattering.particle_);
    Particle tmp_particle2(*particle_);

    vector<CrossSections*> tmp_cross1(scattering.crosssections_);
    vector<CrossSections*> tmp_cross2(crosssections_);

    SetCrossSections(  tmp_cross1 );
    scattering.SetCrossSections(  tmp_cross2 );

    SetParticle( new Particle(tmp_particle1) );
    scattering.SetParticle( new Particle(tmp_particle2) );

    swap(x0_,scattering.x0_);
    swap(do_interpolation_,scattering.do_interpolation_);

    if(scattering.interpolant_ != NULL)
    {
        interpolant_->swap(*scattering.interpolant_);
    }
    else
    {
        interpolant_ = NULL;
    }

    if(scattering.interpolant_diff_ != NULL)
    {
        interpolant_diff_->swap(*scattering.interpolant_diff_) ;
    }
    else
    {
        interpolant_diff_ = NULL;
    }

    integral_->swap(*scattering.integral_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   29.06.2010
*   \author Jan-Hendrik Koehne
*/



//----------------------------------------------------------------------------//


double Scattering::FunctionToIntegral(double energy)
{
    double aux, aux2;
    double result;

    //Same Implentation as in double ProcessCollection::FunctionToIntegral(double energy)
    //It was reimplemented to avoid building a ProcessCollecion Object to calculate or
    //test the scattering class.
    particle_->SetEnergy(energy);
    result  =    0;

    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        aux     =   crosssections_.at(i)->CalculatedEdx();
        result  +=  aux;
    }

    aux = -1/result;
    //End of the reimplementation

    aux2    =   RY*particle_->GetEnergy() / (particle_->GetMomentum() *particle_->GetMomentum());
    aux     *=  aux2*aux2;

    return aux;
}

//----------------------------------------------------------------------------//

double Scattering::CalculateTheta0(double dr, double ei, double ef)
{

    double aux=-1;
    double cutoff=1;
    if(do_interpolation_)
    {
        if(fabs(ei-ef)>fabs(ei)*HALF_PRECISION)
        {
            aux         =   interpolant_->Interpolate(ei);
            double aux2 =   aux - interpolant_->Interpolate(ef);

            if(fabs(aux2)>fabs(aux)*HALF_PRECISION)
            {
                aux =   aux2;
            }
            else
            {
                aux =   interpolant_diff_->Interpolate((ei+ef)/2)*(ef-ei);
            }
        }
        else
        {
            aux =   interpolant_diff_->Interpolate((ei+ef)/2)*(ef-ei);
        }
    }
    else
    {
        aux = integral_->Integrate(ei, ef, boost::bind(&Scattering::FunctionToIntegral, this, _1),4);
    }

    aux =   sqrt(max(aux, 0.0)/x0_) *particle_->GetCharge();
    aux *=  max(1 + 0.038*log(dr/x0_), 0.0);


    return min(aux, cutoff);
}


//----------------------------------------------------------------------------//

double Scattering::FunctionToBuildInterpolant(double energy)
{
        return integral_->Integrate(energy, BIGENERGY, boost::bind(&Scattering::FunctionToIntegral, this, _1),4);
}

void Scattering::EnableInterpolation()
{
    interpolant_ = new Interpolant(NUM2,particle_->GetLow() , BIGENERGY ,boost::bind(&Scattering::FunctionToBuildInterpolant, this, _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );
    interpolant_diff_ = new Interpolant(NUM2,particle_->GetLow() , BIGENERGY ,boost::bind(&Scattering::FunctionToIntegral, this, _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );

    do_interpolation_ = true;
}

void Scattering::DisableInterpolation()
{
    delete interpolant_;
    delete interpolant_diff_;

    interpolant_ = NULL;
    interpolant_diff_ = NULL;

    do_interpolation_ = false;
}

//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Scattering::SetParticle(Particle* particle)
{
    if(particle == NULL || particle == particle_)return;

    for(unsigned int i = 0 ; i< crosssections_.size() ; i++)
    {
        crosssections_.at(i)->SetParticle(particle);
    }

    particle_ = particle;
}

void Scattering::SetCrossSections(std::vector<CrossSections*> crosssections)
{
    crosssections_ = crosssections;
}




