#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Constants.h"
#include <cmath>

Bremsstrahlung::Bremsstrahlung(){

}
//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung &brems)
{
    *this = brems;
}
//----------------------------------------------------------------------------//

Bremsstrahlung& Bremsstrahlung::operator=(const Bremsstrahlung &brems){
    return *this;
}

//----------------------------------------------------------------------------//

void Bremsstrahlung::SetIntegralLimits(){
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableStochasticInerpolation(){
    doStochasticInterpolation_=true;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableContinuousInerpolation(){
    doContinuousInterpolation_=true;
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToContinuousIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToStochasticalIntegral(double variable){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::KelnerKakoulinPetrukhinParametrization(double v, int i)
{

    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow((medium_->GetNucCharge()).at(i), -1./3);

    int step;
    double d, da, dn, Fa, maxV;

    d       =   particle_->GetMass()*particle_->GetMass()
                *v/(2*(particle_->GetEnergy())*(1-v));
    s1      =   (medium_->GetLogConstant()).at(i)*Z3;
    da      =   log(1 + ME/(d*SQRTE*s1));
    Dn      =   1.54*pow((medium_->GetAtomicNum()).at(i), 0.27);
    s1      =   ME*Dn/((particle_->GetMass())*s1);
    dn      =   log(Dn/(1 + d*(Dn*SQRTE - 2)/(particle_->GetMass())));
    maxV    =   ME*(particle_->GetEnergy() - particle_->GetMass())
                /((particle_->GetEnergy())
                  *(particle_->GetEnergy() - particle_->GetMomentum() + ME));

    if(v<maxV)
    {
        Fa  =   log(((particle_->GetMass())/d)/(d*(particle_->GetMass())/(ME*ME) + SQRTE)) -
                log(1 + ME/(d*SQRTE*medium_->GetBPrime().at(i)*(pow(medium_->GetNucCharge().at(i) , -2./3))));
    }
    else
    {
        Fa  =   0;
    }

    if((medium_->GetNucCharge()).at(i)==1)
    {
        step    =   0;
    }

    else
    {
        step    =   1;
    }


    result = ((4./3)*(1-v) + v*v)
            *(log((particle_->GetMass())/d)
              - 0.5 -da - dn + (dn*step + Fa)/(medium_->GetNucCharge().at(i)));

    return result;

}

//----------------------------------------------------------------------------//



double Bremsstrahlung::AndreevBezrukovBugaevParametrization(double v, int i)
{

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;

    Z3 = pow((medium_->GetNucCharge()).at(i), -1./3);

    double aux1, aux2, a1, a2,zeta, qc, qmin, x1, x2, d1,d2, psi1, psi2;

    a1      =   111.7*Z3/ME;
    a2      =   724.2*Z3*Z3/ME;
    qc      =   1.9*MMU*Z3;
    aux     =   2*(particle_->GetMass())/qc;
    aux     *=  aux;
    zeta    =   sqrt(1+aux);
    qmin    =   pow((particle_->GetMass()),2)
                *v/((particle_->GetEnergy())*(1-v));

    x1      =   a1*qmin;
    x2      =   a2*qmin;

    if((medium_->GetNucCharge()).at(i)==1)
    {
        d1  =   0;
        d2  =   0;
    }
    else
    {
        aux1    =   log((particle_->GetMass())/qc);
        aux2    =   (zeta/2)*log((zeta+1)/(zeta-1));
        d1      =   aux1 + aux2;
        d2      =   aux1 + ((3 - pow(zeta , 2))*aux2 + aux)/2;
    }

    aux     =   (particle_->GetMass())*a1;
    aux1    =   log(pow(aux , 2)/(1 + pow(x1 , 2)));
    aux     =   (particle_->GetMass())*a2;
    aux2    =   log(pow(aux , 2)/(1 + pow(x2 , 2)));
    psi1    =   (1+ aux1)/2 + (1 + aux2)/(2*(medium_->GetNucCharge()).at(i));
    psi2    =   (2./3 + aux1)/2 +
                (2./3 + aux2)/(2*(medium_->GetNucCharge()).at(i));

    aux1    =   x1*atan(1/x1);
    aux2    =   x2*atan(1/x2);
    psi1    -=  aux1 + aux2/(medium_->GetNucCharge().at(i));
    aux     =   pow(x1 , 2);
    psi2    +=  2*aux*(1 - aux1 + 3./4*log(aux/(1 + aux)));
    aux     =   pow(x2 , 2);
    psi2    +=  2*aux*(1 - aux2 + 3./4*log(aux/(1 + aux)))
                /(medium_->GetNucCharge().at(i));

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


double Bremsstrahlung::PetrukhinShestakovParametrization(double v, int i)
{

    double Z3       =   0;
    double result   =   0;
    double d, Fd;

    Z3  =   pow((medium_->GetNucCharge()).at(i), -1./3);

    d   =   pow((particle_->GetMass()) , 2)
            * v/(2*(particle_->GetEnergy())*(1-v));

    Fd  =   189*Z3/ME;
    Fd  =   (particle_->GetMass())*Fd/(1 + SQRTE*d*Fd);

    if((medium_->GetNucCharge()).at(i)>10)
    {
        Fd  *=  (2./3)*Z3;
    }

    result  =   ((4./3)*(1-v) + pow(v , 2))*log(Fd);

    return result;

}

//----------------------------------------------------------------------------//



double Bremsstrahlung::CompleteScreeningCase(double v, int i)
{

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Lr, fZ, Lp;

    Z3  =   pow((medium_->GetNucCharge()).at(i) , -1./3);

    aux =   ALPHA*(medium_->GetNucCharge().at(i));
    aux *=  aux;
    fZ  =   aux*(1/(1 + aux) + 0.20206 + aux*(-0.0369 + aux*(0.0083 - 0.002*aux)));

    //check rounding
    switch((int)((medium_->GetNucCharge()).at(i) + 0.5))
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

    result = (((4./3)*(1-v) + pow(v , 2))*
              ((medium_->GetNucCharge()).at(i)*(Lr - fZ) + Lp)
             + (1./9)*(1-v)*((medium_->GetNucCharge()).at(i) + 1))
            /(medium_->GetNucCharge()).at(i);

    return result;

}

//----------------------------------------------------------------------------//

double Bremsstrahlung::ElasticBremsstrahlungCrossSection(double v, int i){

    double aux      =   0;
//    double Z3       =   0;
    double result   =   0;
//    double Dn       =   0;
//    double s1       =   0;

    //Z3  =   pow((medium_->GetNucCharge()).at(i), -1./3);

    switch(parametrization_)
    {
        case 1:
        {
            result  =   KelnerKakoulinPetrukhinParametrization(v, i);
        }break;

        case 2:
        {
            result  =   AndreevBezrukovBugaevParametrization(v, i);
        }break;

        case 3:
        {
            result  =   PetrukhinShestakovParametrization(v, i);
        }break;

        default:
        {
            result  =   CompleteScreeningCase(v, i);
        }break;

    }

    aux =   2*(medium_->GetNucCharge()).at(i)*(ME/particle_->GetMass())*RE;
    aux *=  (ALPHA/v)*aux*result;

//    if(cros->get_lpm())
//    {
//        if(form_!=1)
//        {
//            s1  =   (medium_->GetLogConstant()).at(i)*Z3;
//            Dn  =   1.54*pow((medium_->GetAtomicNum()).at(i) , 0.27);
//            s1  =   ME*Dn/((particle_->GetMass())*s1);
//        }
//        aux *=  lpm(v,s1);
//    }

    double c2   =   pow(particle_->GetCharge() , 2);

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*pow(c2 , 2)*aux;
}
