
#include <cmath>

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/BremsPetrukhinShestakov.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

BremsPetrukhinShestakov::BremsPetrukhinShestakov(const ParticleDef& particle_def,
                                                           const Medium& medium,
                                                           const EnergyCutSettings& cuts,
                                                           Definition param_def)
    : Bremsstrahlung(particle_def, medium, cuts, param_def)
{
}

BremsPetrukhinShestakov::BremsPetrukhinShestakov(const BremsPetrukhinShestakov& brems)
    : Bremsstrahlung(brems)
{
}

// ------------------------------------------------------------------------- //
// AndreevBezrukovBugaev parametrization
// ------------------------------------------------------------------------- //

double BremsPetrukhinShestakov::CalculateParametrization(double energy, double v)
{
    double Z3       =   0;
    double result   =   0;
    double d, Fd;

    Z3  =   pow((current_component_->GetNucCharge()), -1./3);

    d   =   pow((particle_def_.mass) , 2)
            * v/(2*(energy)*(1-v));

    Fd  =   189*Z3/ME;
    Fd  =   (particle_def_.mass)*Fd/(1 + SQRTE*d*Fd);

    if((current_component_->GetNucCharge())>10)
    {
        Fd  *=  (2./3)*Z3;
    }

    result  =   ((4./3)*(1-v) + pow(v , 2))*log(Fd);

    return result;
}
