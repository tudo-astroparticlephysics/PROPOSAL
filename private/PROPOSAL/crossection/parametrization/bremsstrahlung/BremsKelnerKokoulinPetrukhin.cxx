
#include <cmath>

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/BremsKelnerKokoulinPetrukhin.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

BremsKelnerKokoulinPetrukhin::BremsKelnerKokoulinPetrukhin(const ParticleDef& particle_def,
                                                           const Medium& medium,
                                                           const EnergyCutSettings& cuts,
                                                           Definition param_def)
    : Bremsstrahlung(particle_def, medium, cuts, param_def)
{
}

BremsKelnerKokoulinPetrukhin::BremsKelnerKokoulinPetrukhin(const BremsKelnerKokoulinPetrukhin& brems)
    : Bremsstrahlung(brems)
{
}

// ------------------------------------------------------------------------- //
// BremsKelnerKokoulinPetrukhin parametrization
// ------------------------------------------------------------------------- //

double BremsKelnerKokoulinPetrukhin::CalculateParametrization(double energy, double v)
{
    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow(current_component_->GetNucCharge(), -1./3);

    int step;
    double d, da, dn, Fa, maxV;

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(std::max(square_momentum, 0.0));

    d       =   particle_def_.mass*particle_def_.mass
                *v/(2*(energy)*(1-v));
    s1      =   (current_component_->GetLogConstant())*Z3;
    da      =   log(1 + ME/(d*SQRTE*s1));
    Dn      =   1.54*pow((current_component_->GetAtomicNum()), 0.27);
    s1      =   ME*Dn/((particle_def_.mass)*s1);
    dn      =   log(Dn/(1 + d*(Dn*SQRTE - 2)/particle_def_.mass));
    maxV    =   ME*(energy - particle_def_.mass)
                /((energy)
                  *(energy - particle_momentum + ME));

    if(v<maxV)
    {
        Fa  =   log(((particle_def_.mass)/d)/(d*(particle_def_.mass)/(ME*ME) + SQRTE)) -
                log(1 + ME/(d*SQRTE*current_component_->GetBPrime()*(pow(current_component_->GetNucCharge() , -2./3))));
    }
    else
    {
        Fa  =   0;
    }

    if((current_component_->GetNucCharge())==1)
    {
        step    =   0;
    }

    else
    {
        step    =   1;
    }


    result = ((4./3)*(1-v) + v*v)
            *(log((particle_def_.mass)/d)
              - 0.5 -da - dn + (dn*step + Fa)/(current_component_->GetNucCharge()));

    return result;
}
