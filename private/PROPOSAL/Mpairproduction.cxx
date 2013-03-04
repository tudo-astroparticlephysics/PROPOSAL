/*! \file   Mpairproduction.cxx
*   \brief  Source file for the pairproduction routines.
*
*   For more details see the class documentation.
*
*   \date   15.03.2011
*   \author Martin Schmitz
*/

#include "PROPOSAL/Mpairproduction.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/MpairContinuous.h"
#include "PROPOSAL/MpairStochastic.h"




using namespace std;

Mpairproduction::Mpairproduction(Mpairproduction *cros_)
{
    this->particle_ =   cros_->particle_;
    this->medium_   =   cros_->medium_;
    this->cros      =   cros_->cros;
    this->Mpair_    =   cros_;

}

//----------------------------------------------------------------------------//

Mpairproduction::Mpairproduction(CrossSections *cros)
:CrossSections(*cros)
{

    jt_         =   false;
    Mpair_      =   this;
    form_       =   1;
    cm_         =   1;
    continuous_ =   new MpairContinuous(this);
    stochastic_ =   new MpairStochastic(this);
    integral_   =   new Integral(IROMB, IMAXS, IPREC);

}

//----------------------------------------------------------------------------//

double Mpairproduction::mpair(double v, int NumberOfComponent)
{
    if(cm_ <= 0)
    {
        return 0;
    }

    this->v                 =   v;
    this->NumberOfComponent =   NumberOfComponent;

    double aux;
    double Z2;

    if(get_jt())
    {
        setEnergy(NumberOfComponent);

        if(v>=vUp_)
        {
            return max(interpolateJ_.at(NumberOfComponent)->interpolate(particle_->e,log(v/vUp_)/log(vMax_/vUp_)), 0.0);
        }
    }

    if(form_ == 1)
    {
        Z2  =   medium_->get_NucCharge().at(NumberOfComponent)*medium_->get_NucCharge().at(NumberOfComponent);
    }
    else
    {
        Z2  =   medium_->get_NucCharge().at(NumberOfComponent)*(medium_->get_NucCharge().at(NumberOfComponent) + 1);
    }

    rhomax  =   1 - (2*MMU)/(v*particle_->get_energy());
    aux     =   2/(3*PI)*Z2*ALPHA*ALPHA*RM*RM*(1 - v)/v*(integral_->integrateOpened(-rhomax,rhomax,this));
    aux     *=  medium_->get_molDensity()*medium_->get_atomInMolecule().at(NumberOfComponent)*particle_->c*particle_->c;

    return aux;
}

//----------------------------------------------------------------------------//

double Mpairproduction::function(double rho)
{

    if(rho > rhomax || rho < -rhomax)
    {
        return 0;
    }

    double aux;
    double aux2;
    double aux3;
    // rho dependent parts of the equation
    double F;           // Phi
    double X;           // X

    // some constants needed
    double beta;            // beta
    double xi;              // Xi
    double rho2;            // rho^2

    // calculate those constants

    beta    =   v*v/(2*(1-v));
    rho2    =   rho*rho;
    xi      =   v*v*(1-rho*rho)/(4*(1-v));

    // calculate F

    aux     =   ( (2 + rho2)*(1 + beta) + xi*(3 + rho2)  )*log(1 + 1/xi);
    aux2    =   ( (1 + rho2)*(1 + 3/2*beta) - 1/xi*( 1+ 2*beta)*(1 - rho2))*log(1 + xi);
    aux3    =   -1-3*rho2+beta*(1-2*rho2);
    F       =   aux+aux2+aux3;

    // calculate X

    if(form_== 1)
    {
        X   =   1 + calculateU(rho) - calculateU(1 - 2*MMU/(v*particle_->get_energy()));
    }
    else
    {
        X   =   medium_->get_logConstant().at(NumberOfComponent)*MMU/ME;
    }

    return F*log(X);
}

//----------------------------------------------------------------------------//


double Mpairproduction::calculateU(double rho)
{

    double Y;
    double xi;
    double aux, aux2, aux3;

    Y       =   12*sqrt(MMU/particle_->get_energy());
    xi      =   v*v*(1 - rho*rho)/(4*(1 - v));
    aux     =   0.65*pow(medium_->get_atomicNum().at(NumberOfComponent) , -0.27)*medium_->get_logConstant().at(NumberOfComponent)*pow(medium_->get_NucCharge().at(NumberOfComponent) , -1/3)*MMU/ME;
    aux2    =   2*SQRTE*MMU*MMU*medium_->get_logConstant().at(NumberOfComponent)*pow(medium_->get_NucCharge().at(NumberOfComponent) , -1/3)*(1 + xi)*(1 + Y);
    aux3    =   ME*particle_->get_energy()*v*(1-rho*rho);

    return aux/(1 + aux2/aux3);
}

//----------------------------------------------------------------------------//


void Mpairproduction::setEnergy(int i)
{
    cros->set_component(i);

    vMax_   =   1 - MMU/particle_->get_energy();
    vMin_   =   2*MMU/particle_->get_energy();


    if(vMax_<vMin_)
    {
        vMax_   =   vMin_;
    }

    vUp_    =   min(vMax_, medium_->vCut(particle_->e));

    if(vUp_<vMin_)
    {
        vUp_    =   vMin_;
    }
}

//----------------------------------------------------------------------------//

double Mpairproduction::functionInt(double e, double v)
{
    particle_->setEnergy(e);
    setEnergy(get_component());

    if(vUp_==vMax_)
    {
        return 0;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return mpair(v,get_component());
}

//----------------------------------------------------------------------------//

void Mpairproduction::activate_interpolation()
{
    set_jt(false);
    this->interpolateJ_.resize(medium_->get_numCompontents());

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        set_component(i);
        this->interpolateJ_.at(i)   =   new Interpolate(NUM1, e_low, e_hi, NUM1, 0., 1., this, g, false, false, true, g, false, false, false, g, false, false, false);
    }
    set_jt(true);
}



