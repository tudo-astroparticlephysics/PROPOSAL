
#include <cmath>

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/BremsAndreevBezrukovBugaev.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

BremsAndreevBezrukovBugaev::BremsAndreevBezrukovBugaev(const ParticleDef& particle_def,
                                                           const Medium& medium,
                                                           const EnergyCutSettings& cuts,
                                                           Definition param_def)
    : Bremsstrahlung(particle_def, medium, cuts, param_def)
{
}

BremsAndreevBezrukovBugaev::BremsAndreevBezrukovBugaev(const BremsAndreevBezrukovBugaev& brems)
    : Bremsstrahlung(brems)
{
}

// ------------------------------------------------------------------------- //
// AndreevBezrukovBugaev parametrization
// ------------------------------------------------------------------------- //

double BremsAndreevBezrukovBugaev::CalculateParametrization(double energy, double v)
{
    double aux      =   0;
    double Z3       =   0;
    double result   =   0;

    Z3 = pow((current_component_->GetNucCharge()), -1./3);

    double aux1, aux2, a1, a2,zeta, qc, qmin, x1, x2, d1,d2, psi1, psi2;

    a1      =   111.7*Z3/ME;
    a2      =   724.2*Z3*Z3/ME;
    qc      =   1.9*MMU*Z3;
    aux     =   2*(particle_def_.mass)/qc;
    aux     *=  aux;
    zeta    =   sqrt(1+aux);
    qmin    =   pow((particle_def_.mass),2)
                *v/((energy)*(1-v));

    x1      =   a1*qmin;
    x2      =   a2*qmin;

    if((current_component_->GetNucCharge())==1)
    {
        d1  =   0;
        d2  =   0;
    }
    else
    {
        aux1    =   log((particle_def_.mass)/qc);
        aux2    =   (zeta/2)*log((zeta+1)/(zeta-1));
        d1      =   aux1 + aux2;
        d2      =   aux1 + ((3 - pow(zeta , 2))*aux2 + aux)/2;
    }

    aux     =   (particle_def_.mass)*a1;
    aux1    =   log(pow(aux , 2)/(1 + pow(x1 , 2)));
    aux     =   (particle_def_.mass)*a2;
    aux2    =   log(pow(aux , 2)/(1 + pow(x2 , 2)));
    psi1    =   (1+ aux1)/2 + (1 + aux2)/(2*(current_component_->GetNucCharge()));
    psi2    =   (2./3 + aux1)/2 +
                (2./3 + aux2)/(2*(current_component_->GetNucCharge()));

    aux1    =   x1*atan(1/x1);
    aux2    =   x2*atan(1/x2);
    psi1    -=  aux1 + aux2/(current_component_->GetNucCharge());
    aux     =   pow(x1 , 2);
    psi2    +=  2*aux*(1 - aux1 + 3./4*log(aux/(1 + aux)));
    aux     =   pow(x2 , 2);
    psi2    +=  2*aux*(1 - aux2 + 3./4*log(aux/(1 + aux)))
                /(current_component_->GetNucCharge());

    psi1    -=  d1;
    psi2    -=  d2;
    result  =   (2-2*v + pow(v , 2))*psi1 - (2./3)*(1-v)*psi2;

    if(result<0)
    {
        result  =   0;
    }

    return result;
}
