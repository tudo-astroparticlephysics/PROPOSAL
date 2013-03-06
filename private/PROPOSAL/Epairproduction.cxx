/*
 * Epairproduction.cxx
 *
 *  Created on: 29.06.2010
 *      Author: koehne
 */
#include "PROPOSAL/Epairproduction.h"
#include <algorithm>
#include <cmath>
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/EpairContinuous.h"
#include "PROPOSAL/EpairStochastic.h"

using namespace std;


Epairproduction::Epairproduction(Epairproduction *cros_)
{
    init            =   true;
    this->particle_ =   cros_->particle_;
    this->medium_   =   cros_->medium_;
    this->cros      =   cros_->cros;
    epair           =   cros_;

}

//----------------------------------------------------------------------------//


Epairproduction::Epairproduction(CrossSections *cros)
:CrossSections(*cros)
{

init        =   true;
epair       =   this;
continuous_ =   new EpairContinuous(this);
stochastic_ =   new EpairStochastic(this);
integral_   =   new Integral(IROMB, IMAXS, IPREC);
}

//----------------------------------------------------------------------------//

void Epairproduction::setEnergy(int i)
{

    double aux;
    cros->set_component(i);

    vMin    =   4*ME/particle_->e;
    vMax    =   1 - (3./4)*SQRTE*(particle_->m/particle_->e) * pow(medium_->get_NucCharge().at(i) , 1./3);
    aux     =   particle_->m/particle_->e;
    aux     =   1-6*aux*aux;
    vMax    =   min(vMax, aux);
    vMax    =   min(vMax, 1-particle_->m/particle_->e);
    
    if(vMax<vMin)
    {
        vMax    =   vMin;
    }
    
    vUp     =   min(vMax, medium_->vCut(particle_->e));
    
    if(vUp<vMin)
    {
        vUp     =   vMin;
    }
}

//----------------------------------------------------------------------------------------------------//

double Epairproduction::function(double r)
{
    double Fe, Fm, Le, Lm, Ye, Ym, s, b, k, g1, g2;
    double aux, aux1, aux2, r2, Z3;

    r       =   1-r; // only for integral optimization - do not forget to swap integration limits!
    r2      =   r*r;
    Z3      =   pow(medium_->get_NucCharge().at(i) , -1./3);
    aux     =   (particle_->m*v)/(2*ME);
    aux     *=  aux;
    s       =   aux*(1 - r2)/(1 - v);
    b       =   (v*v)/(2*(1 - v));
    Ye      =   (5 - r2 + 4*b*(1 + r2))/(2*(1 + 3*b)*log(3 + 1/s) - r2 - 2*b*(2 - r2));
    Ym      =   (4 + r2 + 3*b*(1 + r2))/((1 + r2)*(1.5 + 2*b)*log(3 + s) + 1 - 1.5*r2);
    aux     =   (1.5*ME)/(particle_->m*Z3);
    aux     *=  aux;
    aux1    =   (1 + s)*(1 + Ye);
    aux2    =   (2*ME*SQRTE*medium_->get_logConstant().at(i)*Z3) / (particle_->e*v*(1 - r2));
    Le      =   log((medium_->get_logConstant().at(i)*Z3*sqrt(aux1)) / (1 + aux2*aux1)) - 0.5*log(1 + aux*aux1);
    Lm      =   log(((particle_->m/(1.5*ME))*medium_->get_logConstant().at(i)*Z3*Z3)/(1 + aux2*(1 + s)*(1 + Ym)));

    if ( Le > 0 )
    {
        if (1/s < HALF_PRECISION)
        {
            Fe = (1.5 - r2/2 + b*(1 + r2))/s*Le;
        }
        else
        {
            Fe = (((2 + r2)*(1 + b) + s*(3 + r2))*log(1 + 1/s) + (1 - r2 - b)/(1 + s) - (3 + r2))*Le;
        }
    }
    else
    {
        Fe = 0;
    }

    if ( Le > 0)
    {
        Fm = (((1 + r2)*(1 + 1.5*b) - (1 + 2*b)*(1 - r2)/s)*log(1 + s) + s*(1 - r2 - b)/(1 + s) + (1 + 2*b)*(1 - r2))*Lm;
    }

    else
    {
        Fm = 0;
    }

    if(medium_->get_NucCharge().at(i)==1)
    {
        g1  =   4.4e-5;
        g2  =   4.8e-5;
    }
    else
    {
        g1  =   1.95e-5;
        g2  =   5.3e-5;
    }

    aux     =   particle_->e/particle_->m;
    k       =   0.058*log(aux/(1 + g2*aux/Z3)) - 0.14;

    if(k<=0)
    {
        k   =   0;
    }

    else
    {
        k   =   (0.073*log(aux/(1 + g1*aux/(Z3*Z3))) - 0.26)/k;
    }

    if(k<0)
    {
        k   =   0;
    }

    aux     =   ALPHA*RE;
    aux     *=  aux/(1.5*PI);
    aux1    =   ME/particle_->m;
    aux1    *=  aux1;

    aux     *=  2*medium_->get_NucCharge().at(i) * (medium_->get_NucCharge().at(i) + k)*((1 - v)/v)*lpm(r2, b, s)*(Fe + aux1*Fm);
    if(aux<0)
    {
        aux =   0;
    }

    return aux;
}

//----------------------------------------------------------------------------//


double Epairproduction::lpm(double r2, double b, double x)
{
    if(cros->get_lpm())
    {
        if(init)
        {
            init        =   false;
            double sum  =   0;

            for(int i=0; i<medium_->get_numCompontents(); i++)
            {
                sum +=  medium_->get_NucCharge().at(i)*medium_->get_NucCharge().at(i)*log(3.25*medium_->get_logConstant().at(i)*pow(medium_->get_NucCharge().at(i), -1./3));
            }

            eLpm    =   particle_->m/(ME*RE);
            eLpm    *=  (eLpm*eLpm)*ALPHA*particle_->m/(2*PI*medium_->get_molDensity()*particle_->c*particle_->c*sum);
        }

        double A, B, C, D, E, s;
        double s2, s36, s6, d1, d2, atan_, log1, log2;

        s       =   sqrt(eLpm/(particle_->e*v*(1 - r2)))/4;
        s6      =   6*s;
        atan_   =   s6*(x + 1);

        if(atan_>1/COMPUTER_PRECISION)
        {
            return 1;
        }

        s2      =   s*s;
        s36     =   36*s2;
        d1      =   s6/(s6 + 1);
        d2      =   s36/(s36 + 1);
        atan_   =   atan(atan_) - PI/2;
        log1    =   log((s36*(1 + x)*(1 + x) + 1)/(s36*x*x));
        log2    =   log((s6*(1 + x) + 1)/(s6*x));
        A       =   0.5*d2*(1 + 2*d2*x)*log1 - d2 + 6*d2*s*(1 + ((s36 - 1)/(s36 + 1))*x)*atan_;
        B       =   d1*(1 + d1*x)*log2 - d1;
        C       =   -d2*d2*x*log1 + d2 - (d2*d2*(s36 - 1)/(6*s))*x*atan_;
        D       =   d1-d1*d1*x*log2;
        E       =   -s6*atan_;

        return ((1 + b)*(A + (1 + r2)*B) + b*(C + (1 + r2)*D) + (1 - r2)*E)/(((2 + r2)*(1 +b ) + x*(3 + r2))*log(1 +1 /x) + (1 - r2 - b)/(1 + x) - (3 + r2));
    }

    else
    {
        return 1;
    }
}

//----------------------------------------------------------------------------//



double Epairproduction::ePair(double v, int i)
{

    if(jt_)
    {
        setEnergy(i);

        if(v>=vUp)
        {
            return max(interpolateJ_[i].interpolate(particle_->e, log(v/vUp)/log(vMax/vUp)), 0.0);
        }
    }

    double rMax, aux, aux2;

    this->i =   i;
    this->v =   v;
    aux     =   1 - (4*ME)/(particle_->e*v);
    aux2    =   1 - (6*particle_->m*particle_->m)/(particle_->e*particle_->e*(1 - v));

    if(aux>0 && aux2>0)
    {
        rMax    =   sqrt(aux)*aux2;
    }
    else
    {
        rMax    =   0;
    }

    aux =   max(1 - rMax , COMPUTER_PRECISION);

    return (1 + epairerror)*medium_->get_molDensity()*medium_->get_atomInMolecule().at(i)*particle_->c*particle_->c*(integral_->integrateOpened(1 - rMax, aux, this) + integral_->integrateWithLog(aux, 1, this));

}

//----------------------------------------------------------------------------//

double Epairproduction::functionInt(double e, double v)
{
    particle_->setEnergy(e);
    setEnergy(cros->get_component());

    if(vUp==vMax)
    {
    return 0;
    }

    v   =   vUp*exp(v*log(vMax/vUp));

    return ePair(v, cros->get_component());
}

void Epairproduction::SetRandomNumberGenerator(boost::function<double ()> &f)
{
	MathModel::SetRandomNumberGenerator(f);
	get_Continuous()->SetRandomNumberGenerator(f);
	get_Stochastic()->SetRandomNumberGenerator(f);
}


