
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/Constants.h"

// #define BREMSSTRAHLUNG_IMPL(param)                                                                                     \
//     Brems##param::Brems##param(                                                                                        \
//         const ParticleDef& particle_def, const Medium& medium, const EnergyCutSettings& cuts, Definition param_def)    \
//         : Bremsstrahlung(particle_def, medium, cuts, param_def)                                                        \
//     {                                                                                                                  \
//     }                                                                                                                  \
//                                                                                                                        \
//     Brems##param::Brems##param(const Brems##param& brems)                                                              \
//         : Bremsstrahlung(brems)                                                                                        \
//     {                                                                                                                  \
//     }                                                                                                                  \
//                                                                                                                        \
//     Brems##param::~Brems##param() {}

using namespace PROPOSAL;

/******************************************************************************
*                                   HardBB                                    *
******************************************************************************/

std::vector<const double> HardBB::x = boost::assign::list_of(3)(4)(5)(6)(7)(8)(9);

HardBB::HardBB(const ParticleDef& particle_def)
    : particle_def_(particle_def)
    , interpolant_()
{
    const HardBBTables::VecType* y = particle_def_.hardbb_table;

    if (y != NULL)
    {
        for(unsigned int i=0; i < y->size(); i++)
        {
            interpolant_.push_back(new Interpolant(x, y->at(i), 4, false, false));
        }
    }
}

HardBB::HardBB(const HardBB& hardBB)
    : particle_def_(hardBB.particle_def_)
    , interpolant_()
{
    for(unsigned int i = 0; i < hardBB.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*hardBB.interpolant_[i]);
    }
}

HardBB::~HardBB()
{
    for(std::vector<Interpolant*>::iterator it = interpolant_.begin(); it != interpolant_.end(); ++it)
    {
        delete *it;
    }
}

// ------------------------------------------------------------------------- //
double HardBB::CalculateHardBB(double energy, double v)
{
    if (energy < 1.0e5 || v < 1.0e-7)
    {
        return 0;
    }

    double aux, sum, lov, loe;

    sum = 0;
    aux = 1;
    lov = log(v) / LOG10;
    loe = log(energy) / LOG10 - 3;

    for (unsigned int i = 0; i < interpolant_.size(); i++)
    {
        if (i > 0)
        {
            aux *= lov;
        }

        sum += aux * interpolant_[i]->InterpolateArray(loe);
    }
    return sum / v;
}

/******************************************************************************
*                                ShadowEffect                                *
******************************************************************************/

// ------------------------------------------------------------------------- //
double ShadowDutta::CalculateShadowEffect(const Components::Component& component, double x, double nu)
{
    (void) nu;

    if (component.GetNucCharge() == 1)
        return 1;

    double G;

    if(x<0.0014)
    {
        G   =   pow(component.GetAtomicNum(), -0.1);
    }
    else if(x<0.04)
    {
        G   =   pow(component.GetAtomicNum(), 0.069*log(x)/LOG10+0.097);
    }
    else
    {
        G   =   1;
    }

    return G;
}

// ------------------------------------------------------------------------- //
double ShadowEffect::CalculateShadowEffect(const Components::Component& component, double x, double nu)
{
    if(component.GetNucCharge()==1) return 1;

    double G;

    if(x>0.3)
    {
        const double Mb =   0.437;
        const double la =   0.5;
        const double x2 =   0.278;

        double mb, Aosc, mu, au, ac;

        mb      =   Mb*component.GetMN();
        au      =   1/(1 - x);
        ac      =   1/(1 - x2);
        mu      =   MPI/component.GetAverageNucleonWeight();
        Aosc    =   (1 - la*x)*((au - ac)-mu*(au*au - ac*ac));
        G       =   1 - mb*Aosc;
    }
    else
    {
        const double M1 =   0.129;
        const double M2 =   0.456;
        const double M3 =   0.553;

        double m1, m2, m3, x0, sgn, tmp;

        m1  =   M1*component.GetMN();
        m2  =   M2*component.GetMN();
        m3  =   M3*component.GetMN();
        nu  *=  1.e-3;
        sgn =   112.2*(0.609*pow(nu, 0.0988) + 1.037*pow(nu, -0.5944));

        // Bezrukav Bugaev shadow
        tmp =   0.00282*pow(component.GetAtomicNum(), 1./3)*sgn;
        G   =   (3/tmp)*(0.5 + ((1 + tmp)*exp(-tmp) - 1)/(tmp*tmp));

        G   =   0.75*G + 0.25;
        x0  =   pow(G/(1+m2), 1/m1);

        if(x>=x0)
        {
            G   =   pow(x, m1)*(1+m2)*(1-m3*x);
        }
    }

    return G;
}

/******************************************************************************
*                               Photonuclear                                  *
******************************************************************************/


// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Photonuclear::Photonuclear(const ParticleDef& particle_def,
                           const Medium& medium,
                           const EnergyCutSettings& cuts,
                           Definition param_def)
    : Parametrization(particle_def, medium, cuts, param_def)
{
}

Photonuclear::Photonuclear(const Photonuclear& brems)
    : Parametrization(brems)
{
}

Photonuclear::~Photonuclear()
{
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double Photonuclear::DifferentialCrossSection(double energy, double v)
{
    double aux      =   0;
    double result   =   0;

    result = CalculateParametrization(energy, v);

    aux =   2*(components_[component_index_]->GetNucCharge())*(ME/particle_def_.mass)*RE;
    aux *=  (ALPHA/v)*aux*result;

    if(param_def_.lpm_effect_enabled)
    {
        aux *=  lpm(energy, v);
    }

    double c2   =   pow(particle_def_.charge , 2);

    return medium_->GetMolDensity()*components_[component_index_]->GetAtomInMolecule()*c2*c2*aux;
}

// ------------------------------------------------------------------------- //
double Photonuclear::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * DifferentialCrossSection(energy, variable);
}

//----------------------------------------------------------------------------//
double Photonuclear::FunctionToDNdxIntegral(double energy, double variable)
{
    return param_def_.multiplier * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Photonuclear::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 0.;
    limits.vMax =
        1 - (3. / 4) * SQRTE * (particle_def_.mass / energy) * pow((components_[component_index_]->GetNucCharge()), 1. / 3);

    if (limits.vMax < 0)
    {
        limits.vMax = 0;
    }

    if (lorenz_)
    {
        limits.vMax = std::min(limits.vMax, lorenz_cut_ / (energy));
    }

    limits.vMax = std::min(limits.vMax, (1 - (particle_def_.mass / energy)));
    limits.vUp  = std::min(limits.vMax, cut_settings_.GetCut(energy));

    return limits;
}

/******************************************************************************
*                          Specifc Parametrizations                           *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

// BREMSSTRAHLUNG_IMPL(PetrukhinShestakov)
// BREMSSTRAHLUNG_IMPL(KelnerKokoulinPetrukhin)
// BREMSSTRAHLUNG_IMPL(CompleteScreening)
// BREMSSTRAHLUNG_IMPL(AndreevBezrukovBugaev)

// ------------------------------------------------------------------------- //
// PetrukhinShestakov parametrization
// ------------------------------------------------------------------------- //

double BremsPetrukhinShestakov::CalculateParametrization(double energy, double v)
{
    double Z3       =   0;
    double result   =   0;
    double d, Fd;

    Z3  =   pow((components_[component_index_]->GetNucCharge()), -1./3);

    d   =   pow((particle_def_.mass) , 2)
            * v/(2*(energy)*(1-v));

    Fd  =   189*Z3/ME;
    Fd  =   (particle_def_.mass)*Fd/(1 + SQRTE*d*Fd);

    if((components_[component_index_]->GetNucCharge())>10)
    {
        Fd  *=  (2./3)*Z3;
    }

    result  =   ((4./3)*(1-v) + pow(v , 2))*log(Fd);

    return result;
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

    Z3  =   pow(components_[component_index_]->GetNucCharge(), -1./3);

    int step;
    double d, da, dn, Fa, maxV;

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(std::max(square_momentum, 0.0));

    d       =   particle_def_.mass*particle_def_.mass
                *v/(2*(energy)*(1-v));
    s1      =   (components_[component_index_]->GetLogConstant())*Z3;
    da      =   log(1 + ME/(d*SQRTE*s1));
    Dn      =   1.54*pow((components_[component_index_]->GetAtomicNum()), 0.27);
    dn      =   log(Dn/(1 + d*(Dn*SQRTE - 2)/particle_def_.mass));
    maxV    =   ME*(energy - particle_def_.mass)
                /((energy)
                  *(energy - particle_momentum + ME));

    if(v<maxV)
    {
        Fa  =   log(((particle_def_.mass)/d)/(d*(particle_def_.mass)/(ME*ME) + SQRTE)) -
                log(1 + ME/(d*SQRTE*components_[component_index_]->GetBPrime()*(pow(components_[component_index_]->GetNucCharge() , -2./3))));
    }
    else
    {
        Fa  =   0;
    }

    if((components_[component_index_]->GetNucCharge())==1)
    {
        step    =   0;
    }

    else
    {
        step    =   1;
    }


    result = ((4./3)*(1-v) + v*v)
            *(log((particle_def_.mass)/d)
              - 0.5 -da - dn + (dn*step + Fa)/(components_[component_index_]->GetNucCharge()));

    return result;
}

// ------------------------------------------------------------------------- //
// CompleteScreening parametrization
// ------------------------------------------------------------------------- //

double BremsCompleteScreening::CalculateParametrization(double energy, double v)
{
    (void) energy;

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Lr, fZ, Lp;

    Z3  =   pow((components_[component_index_]->GetNucCharge()) , -1./3);

    aux =   ALPHA*(components_[component_index_]->GetNucCharge());
    aux *=  aux;
    fZ  =   aux*(1/(1 + aux) + 0.20206 + aux*(-0.0369 + aux*(0.0083 - 0.002*aux)));

    //check rounding
    switch((int)((components_[component_index_]->GetNucCharge()) + 0.5))
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
              ((components_[component_index_]->GetNucCharge())*(Lr - fZ) + Lp)
             + (1./9)*(1-v)*((components_[component_index_]->GetNucCharge()) + 1))
            /(components_[component_index_]->GetNucCharge());

    return result;
}

// ------------------------------------------------------------------------- //
// AndreevBezrukovBugaev parametrization
// ------------------------------------------------------------------------- //

double BremsAndreevBezrukovBugaev::CalculateParametrization(double energy, double v)
{
    double aux      =   0;
    double Z3       =   0;
    double result   =   0;

    Z3 = pow((components_[component_index_]->GetNucCharge()), -1./3);

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

    if((components_[component_index_]->GetNucCharge())==1)
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
    psi1    =   (1+ aux1)/2 + (1 + aux2)/(2*(components_[component_index_]->GetNucCharge()));
    psi2    =   (2./3 + aux1)/2 +
                (2./3 + aux2)/(2*(components_[component_index_]->GetNucCharge()));

    aux1    =   x1*atan(1/x1);
    aux2    =   x2*atan(1/x2);
    psi1    -=  aux1 + aux2/(components_[component_index_]->GetNucCharge());
    aux     =   pow(x1 , 2);
    psi2    +=  2*aux*(1 - aux1 + 3./4*log(aux/(1 + aux)));
    aux     =   pow(x2 , 2);
    psi2    +=  2*aux*(1 - aux2 + 3./4*log(aux/(1 + aux)))
                /(components_[component_index_]->GetNucCharge());

    psi1    -=  d1;
    psi2    -=  d2;
    result  =   (2-2*v + pow(v , 2))*psi1 - (2./3)*(1-v)*psi2;

    if(result<0)
    {
        result  =   0;
    }

    return result;
}

// #undef BREMSSTRAHLUNG_IMPL
