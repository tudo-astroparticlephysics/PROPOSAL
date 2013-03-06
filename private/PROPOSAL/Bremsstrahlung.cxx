/*
 * Bremsstrahlung.cxx
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/BremsStochastic.h"
#include "PROPOSAL/BremsContinuous.h"
#include "PROPOSAL/Medium.h"
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------//

// constructors

Bremsstrahlung::Bremsstrahlung()
:vMax_(0)
,vUp_ (0)
,vMin_(0)
,form_(1)
,lorenzCut_(1.e6)
{
    lorenz_ =   false;
    init_   =   true;
}

//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(Bremsstrahlung *cros)
:vMax_(0)
,vUp_ (0)
,vMin_(0)
,form_(1)
,lorenzCut_(1.e6)
{
    lorenz_     =   false;
    init_       =   true;
    particle_   =   cros->particle_;
    medium_     =   cros->medium_;
    this->cros  =   cros->cros;
    brems_      =   cros;
}

//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(CrossSections *cros)
:CrossSections(*cros)
,vMax_(0)
,vUp_ (0)
,vMin_(0)
,form_(1)
,lorenzCut_(1.e6)
{
    lorenz_     =   false;
    init_       =   true;
    brems_      =   this;
    continuous_ =   new BremsContinuous(this);
    stochastic_ =   new BremsStochastic(this);
    integral_   =   new Integral(IROMB, IMAXS, IPREC);

}

//----------------------------------------------------------------------------//

// destructors


Bremsstrahlung::~Bremsstrahlung()
{

}

//----------------------------------------------------------------------------//

// Memberfunctions


void Bremsstrahlung::setEnergy(int i)
{

    cros->set_component(i);
    vMax_   =   1 - (3./4)*SQRTE*(particle_->get_mass()/particle_->get_energy())*pow((medium_->get_NucCharge().at(i)) , 1./3);

    if(vMax_<0)
    {
        vMax_   =   0;
    }

    if(lorenz_)
    {
        vMax_   =   min(vMax_, lorenzCut_/(particle_->get_energy()));
    }

    vMax_   =   min(vMax_, (1-(particle_->get_mass()/particle_->get_energy())));
    vUp_    =   min(vMax_, (medium_->vCut(particle_->get_energy())));
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::Sel(double v, int i)
{
    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow((medium_->get_NucCharge()).at(i), -1./3);

    switch(form_)
    {
        case 1:
        {
            result  =   Kelner_Kakoulin_Petrukhin_parametrization(v, i);
        }break;

        case 2:
        {
            result  =   Andreev_Bezrukov_Bugaev_parametrization(v, i);
        }break;

        case 3:
        {
            result  =   Petrukhin_Shestakov_parametrization(v, i);
        }break;

        default:
        {
            result  =   complete_screening_case(v, i);
        }break;

    }

    aux =   2*(medium_->get_NucCharge()).at(i)*(ME/particle_->get_mass())*RE;
    aux *=  (ALPHA/v)*aux*result;

    if(cros->get_lpm())
    {
        if(form_!=1)
        {
            s1  =   (medium_->get_logConstant()).at(i)*Z3;
            Dn  =   1.54*pow((medium_->get_atomicNum()).at(i) , 0.27);
            s1  =   ME*Dn/((particle_->get_mass())*s1);
        }
        aux *=  lpm(v,s1);
    }

    double c2   =   pow(particle_->get_charge() , 2);

    return (1 + bremserror)*(medium_->get_molDensity())*medium_->get_atomInMolecule().at(i)*pow(c2 , 2)*aux;

}
	
//----------------------------------------------------------------------------//

void Bremsstrahlung::setLpm()
{

    if(init_){

        bool lpmSave    =   cros->get_lpm();
        init_           =   false;

        cros->set_lpm(false);

        double sum      =   0;
        double e        =   particle_->get_energy();

        particle_->setEnergy(BIGENERGY);

        for(int i=0; i < medium_->get_numCompontents(); i++)
        {
            setEnergy(i);
            sum +=  integral_->integrateOpened(0, vUp_, continuous_) + integral_->integrateWithLog(vUp_, vMax_, continuous_);
        }

        double c2   =   pow(particle_->get_charge() , 2);
        xo_         =   pow(c2,2)/sum;
        eLpm_       =   ALPHA*(particle_->get_mass());
        eLpm_       *=  eLpm_/(4*PI*ME*RE*sum);

        particle_->setEnergy(e);
        setEnergy(0);
        cros->set_lpm(lpmSave);

    }
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::lpm(double v, double s1)
{

    if(cros->get_lpm())
    {
        setLpm();
        double G, fi, xi, sp, h, s, s2, s3, ps, Gamma;

        const double fi1    =   1.54954;
        const double G1     =   0.710390;
        const double G2     =   0.904912;
        s1                  *=  s1*SQRT2;
        sp                  =   sqrt(eLpm_*v/(8*(particle_->get_energy())*(1-v)));
        h                   =   log(sp)/log(s1);

        if(sp < s1)
        {
            xi  =   2;
        }
        else if(sp < 1)
        {
            xi  =   1 + h - 0.08*(1 - h)*(1 - (1-h)*(1 - h))/log(s1);
        }
        else
        {
            xi  =   1;
        }

        s       =   sp/sqrt(xi);
        Gamma   =   RE*ME/(ALPHA*(particle_->get_mass())*v);
        Gamma   =   1 +4*PI*(medium_->get_molDensity())*(medium_->get_sumCharge())*RE*pow(Gamma,2);
        s       *=  Gamma;
        s2      =   pow(s,2);
        s3      =   pow(s,3);

        if(s < fi1)
        {
            fi  =   1-exp(-6*s*(1 + (3-PI)*s) + s3/(0.623 + 0.796*s + 0.658*s2));
        }
        else
        {
            fi  =   1 - 0.012/pow(s2 , 2);
        }

        if(s < G1)
        {
            ps  =   1 - exp(-4*s - 8*s2/(1 + 3.936*s + 4.97*s2 - 0.05*s3 + 7.50*pow(s2 , 2)));
            G   =   3*ps - 2*fi;
        }
        else if (s < G2)
        {
            G   =   36*s2/(36*s2 + 1);
        }
        else
        {
            G   =   1 - 0.022/pow(s2 , 2);
        }

        return ((xi/3)*((v*v)*G/(Gamma*Gamma) + 2*(1 + (1-v)*(1-v))*fi/Gamma))/((4./3)*(1-v) + v*v);

    }

    else
    {
        return 1;
    }

}

//----------------------------------------------------------------------------//


double Bremsstrahlung::Kelner_Kakoulin_Petrukhin_parametrization(double v, int i)
{

    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow((medium_->get_NucCharge()).at(i), -1./3);

    int step;
    double d, da, dn, Fa, maxV;

    d       =   particle_->get_mass()*particle_->get_mass()*v/(2*(particle_->get_energy())*(1-v));
    s1      =   (medium_->get_logConstant()).at(i)*Z3;
    da      =   log(1 + ME/(d*SQRTE*s1));
    Dn      =   1.54*pow((medium_->get_atomicNum()).at(i), 0.27);
    s1      =   ME*Dn/((particle_->get_mass())*s1);
    dn      =   log(Dn/(1 + d*(Dn*SQRTE - 2)/(particle_->get_mass())));
    maxV    =   ME*(particle_->get_energy() - particle_->get_mass())/((particle_->get_energy())*(particle_->get_energy() - particle_->get_momentum() + ME));

    if(v<maxV)
    {
        Fa  =   log(((particle_->get_mass())/d)/(d*(particle_->get_mass())/(ME*ME) + SQRTE)) -
        log(1 + ME/(d*SQRTE*medium_->get_bPrime().at(i)*(pow(medium_->get_NucCharge().at(i) , -2./3))));
    }
    else
    {
        Fa  =   0;
    }

    if((medium_->get_NucCharge()).at(i)==1)
    {
        step    =   0;
    }

    else
    {
        step    =   1;
    }


    result = ((4./3)*(1-v) + v*v)*(log((particle_->get_mass())/d) - 0.5 -da - dn + (dn*step + Fa)/(medium_->get_NucCharge().at(i)));

    return result;

}

//----------------------------------------------------------------------------//



double Bremsstrahlung::Andreev_Bezrukov_Bugaev_parametrization(double v, int i)
{

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;

    Z3 = pow((medium_->get_NucCharge()).at(i), -1./3);

    double aux1, aux2, a1, a2,zeta, qc, qmin, x1, x2, d1,d2, psi1, psi2;

    a1      =   111.7*Z3/ME;
    a2      =   724.2*Z3*Z3/ME;
    qc      =   1.9*MMU*Z3;
    aux     =   2*(particle_->get_mass())/qc;
    aux     *=  aux;
    zeta    =   sqrt(1+aux);
    qmin    =   pow((particle_->get_mass()),2)*v/((particle_->get_energy())*(1-v));
    x1      =   a1*qmin;
    x2      =   a2*qmin;

    if((medium_->get_NucCharge()).at(i)==1)
    {
        d1  =   0;
        d2  =   0;
    }
    else
    {
        aux1    =   log((particle_->get_mass())/qc);
        aux2    =   (zeta/2)*log((zeta+1)/(zeta-1));
        d1      =   aux1 + aux2;
        d2      =   aux1 + ((3 - pow(zeta , 2))*aux2 + aux)/2;
    }

    aux     =   (particle_->get_mass())*a1;
    aux1    =   log(pow(aux , 2)/(1 + pow(x1 , 2)));
    aux     =   (particle_->get_mass())*a2;
    aux2    =   log(pow(aux , 2)/(1 + pow(x2 , 2)));
    psi1    =   (1+ aux1)/2 + (1 + aux2)/(2*(medium_->get_NucCharge()).at(i));
    psi2    =   (2./3 + aux1)/2 + (2./3 + aux2)/(2*(medium_->get_NucCharge()).at(i));
    aux1    =   x1*atan(1/x1);
    aux2    =   x2*atan(1/x2);
    psi1    -=  aux1 + aux2/(medium_->get_NucCharge().at(i));
    aux     =   pow(x1 , 2);
    psi2    +=  2*aux*(1 - aux1 + 3./4*log(aux/(1 + aux)));
    aux     =   pow(x2 , 2);
    psi2    +=  2*aux*(1 - aux2 + 3./4*log(aux/(1 + aux)))/(medium_->get_NucCharge().at(i));
    psi1    -=  d1;
    psi2    -=  d2;
    result  =   (2-2*v + pow(v , 2))*psi1 - (2./3)*(1-v)*psi2;

    if(result<0)
    {
        result  =   0;
    }

    return result;

}


//----------------------------------------------------------------------------//


double Bremsstrahlung::Petrukhin_Shestakov_parametrization(double v, int i)
{

    double Z3       =   0;
    double result   =   0;
    double d, Fd;

    Z3  =   pow((medium_->get_NucCharge()).at(i), -1./3);

    d   =   pow((particle_->get_mass()) , 2) * v/(2*(particle_->get_energy())*(1-v));
    Fd  =   189*Z3/ME;
    Fd  =   (particle_->get_mass())*Fd/(1 + SQRTE*d*Fd);

    if((medium_->get_NucCharge()).at(i)>10)
    {
        Fd  *=  (2./3)*Z3;
    }

    result  =   ((4./3)*(1-v) + pow(v , 2))*log(Fd);

    return result;

}
        
//----------------------------------------------------------------------------//


        
double Bremsstrahlung::complete_screening_case(double v, int i)
{

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Lr, fZ, Lp;

    Z3  =   pow((medium_->get_NucCharge()).at(i) , -1./3);

    aux =   ALPHA*(medium_->get_NucCharge().at(i));
    aux *=  aux;
    fZ  =   aux*(1/(1 + aux) + 0.20206 + aux*(-0.0369 + aux*(0.0083 - 0.002*aux)));

    //check rounding
    switch((int)((medium_->get_NucCharge()).at(i) + 0.5))
    {

        case 1:
        {
            Lr  =   5.31;
            Lp  =   6.144;
        }break;

        case 2:
        {
            Lr  =   4.79;
            Lp  =   5.621;
        }break;

        case 3:
        {
            Lr  =   4.74;
            Lp  =   5.805;
        }break;

        case 4:
        {
            Lr  =   4.71;
            Lp  =   5.924;
        }break;

        default:
        {
            Lr  =   log(184.15*Z3);
            Lp  =   log (1194*pow(Z3 , 2));
        }break;

    }

    result = (((4./3)*(1-v) + pow(v , 2))*((medium_->get_NucCharge()).at(i)*(Lr - fZ) + Lp) + (1./9)*(1-v)*((medium_->get_NucCharge()).at(i) + 1))/(medium_->get_NucCharge()).at(i);

    return result;

}

//----------------------------------------------------------------------------//

// Getter

BremsContinuous* Bremsstrahlung::get_Continuous()
{
    return continuous_;
}

BremsStochastic* Bremsstrahlung::get_Stochastic()
{
    return stochastic_;
}

//----------------------------------------------------------------------------//

// Setter

void Bremsstrahlung::set_continous(BremsContinuous *continous)
{
    continuous_ =   continous;
}

void Bremsstrahlung::set_stochastic(BremsStochastic *stochastic)
{
    stochastic_ =   stochastic;
}

void Bremsstrahlung::set_vMax(double vMax)
{
    vMax_   =   vMax;
}

void Bremsstrahlung::set_vUp(double vUp)
{
    vUp_    =   vUp;
}

void Bremsstrahlung::set_vMin(double vMin)
{
    vMin_   =   vMin;
}

void Bremsstrahlung::set_integral(Integral *integral)
{
    integral_   =   integral;
}

void Bremsstrahlung::set_form(int form)
{
    form_   =   form;
}

void Bremsstrahlung::set_init(bool init)
{
    init_   =   init;
}

void Bremsstrahlung::set_eLpm(double eLpm)
{
    eLpm_   =   eLpm;
}

void Bremsstrahlung::set_xo(double xo)
{
    xo_ =   xo;
}

void Bremsstrahlung::set_lorenz(bool lorenz)
{
    lorenz_ =   lorenz;
}

void Bremsstrahlung::set_lorenzCut(double lorenzCut)
{
    lorenzCut_  =   lorenzCut;
}

void Bremsstrahlung::SetRandomNumberGenerator(boost::function<double ()> &f)
{
	MathModel::SetRandomNumberGenerator(f);
	get_Continuous()->SetRandomNumberGenerator(f);
	get_Stochastic()->SetRandomNumberGenerator(f);
}


