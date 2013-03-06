/*! \file   Ionizationloss.cxx
*   \brief  Source file for the ionizationloss routines.
*
*   For more details see the class documentation.
*
*   \date   22.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/IonizContinuous.h"
#include "PROPOSAL/IonizStochastic.h"
#include "PROPOSAL/Ionizationloss.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Medium.h"
#include <cmath>
#include <algorithm>


using namespace std;


Ionizationloss::Ionizationloss(Ionizationloss *cros){

    this->particle_ =   cros->particle_;
    this->medium_   =   cros->medium_;
    this->cros      =   cros->cros;
    this->ioniz     =   cros;
}

//----------------------------------------------------------------------------//



Ionizationloss::Ionizationloss(CrossSections *cros)
:CrossSections(*cros)
{

    this->particle_ =   cros->particle_;
    this->medium_   =   cros->medium_;
    this->cros      =   cros;
    ioniz           =   this;
    cont            =   new IonizContinuous(this);
    stochastic_     =   new IonizStochastic(this);

}

//----------------------------------------------------------------------------//


void Ionizationloss::setEnergy()
{

    double aux;

    beta    =   particle_->p/particle_->e;
    gamma   =   particle_->e/particle_->m;
    vMin    =   (1.e-6*medium_->get_I())/particle_->e;
    aux     =   ME/particle_->m;
    vMax    =   2*ME*(gamma*gamma-1)/((1 + 2*gamma*aux + aux*aux)*particle_->e);
    vMax    =   min(vMax, 1 - particle_->m/particle_->e);

    if(vMax<vMin)
    {
        vMax    =   vMin;
    }

    vUp =   min(vMax, medium_->vCut(particle_->e));

    if(vUp<vMin)
    {
        vUp=vMin;
    }

    if(ioniz!=this)
    {
        ioniz->beta     =   beta;
        ioniz->gamma    =   gamma;
        ioniz->vMax     =   vMax;
    }
}

//----------------------------------------------------------------------------//


double Ionizationloss::d2Ndvdx(double v)
{
    double result, aux, aux2;

    aux     =   beta*beta;
    aux2    =   v/(1 + 1/gamma);
    aux2    *=  0.5*aux2;
    result  =   1 - aux*(v/vMax) + aux2;
    result  *=  IONK*particle_->c*particle_->c*medium_->get_ZA()/(2*aux*particle_->e*v*v);

    return (1 + bremserror)*medium_->get_massDensity()*result;
}

//----------------------------------------------------------------------------//



double Ionizationloss::inelCorrection(double v)
{
    double result, a, b, c;

    a       =   log(1 + 2*v*particle_->e/ME);
    b       =   log((1 - v/vMax)/(1 - v));
    c       =   log((2*gamma*(1 - v)*ME)/(particle_->m*v));
    result  =   a*(2*b + c) - b*b;

    return (ALPHA/(2*PI))*result;
}
//----------------------------------------------------------------------------//

// Getter

IonizContinuous* Ionizationloss::get_Continuous()
{
    return cont;
}

IonizStochastic* Ionizationloss::get_Stochastic()
{
    return stochastic_;
}

void Ionizationloss::SetRandomNumberGenerator(boost::function<double ()> &f)
{
	MathModel::SetRandomNumberGenerator(f);
	get_Continuous()->SetRandomNumberGenerator(f);
	get_Stochastic()->SetRandomNumberGenerator(f);
}
