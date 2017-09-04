
#include <cmath>

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/BremsCompleteScreening.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

BremsCompleteScreening::BremsCompleteScreening(const ParticleDef& particle_def,
                                                           const Medium& medium,
                                                           const EnergyCutSettings& cuts,
                                                           Definition param_def)
    : Bremsstrahlung(particle_def, medium, cuts, param_def)
{
}

BremsCompleteScreening::BremsCompleteScreening(const BremsCompleteScreening& brems)
    : Bremsstrahlung(brems)
{
}

// ------------------------------------------------------------------------- //
// AndreevBezrukovBugaev parametrization
// ------------------------------------------------------------------------- //

double BremsCompleteScreening::CalculateParametrization(double energy, double v)
{
    (void) energy;

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Lr, fZ, Lp;

    Z3  =   pow((current_component_->GetNucCharge()) , -1./3);

    aux =   ALPHA*(current_component_->GetNucCharge());
    aux *=  aux;
    fZ  =   aux*(1/(1 + aux) + 0.20206 + aux*(-0.0369 + aux*(0.0083 - 0.002*aux)));

    //check rounding
    switch((int)((current_component_->GetNucCharge()) + 0.5))
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
              ((current_component_->GetNucCharge())*(Lr - fZ) + Lp)
             + (1./9)*(1-v)*((current_component_->GetNucCharge()) + 1))
            /(current_component_->GetNucCharge());

    return result;
}
