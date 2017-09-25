
#include <boost/bind.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/Constants.h"

#define BREMSSTRAHLUNG_IMPL(param)                                                                                     \
    Brems##param::Brems##param(const ParticleDef& particle_def,                                                        \
                               const Medium& medium,                                                                   \
                               const EnergyCutSettings& cuts,                                                          \
                               double multiplier,                                                                      \
                               bool lpm)                                                                               \
        : Bremsstrahlung(particle_def, medium, cuts, multiplier, lpm)                                                  \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Brems##param::Brems##param(const Brems##param& brems)                                                              \
        : Bremsstrahlung(brems)                                                                                        \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Brems##param::~Brems##param() {}                                                                                   \
                                                                                                                       \
    const std::string Brems##param::name_ = "Brems" #param;

using namespace PROPOSAL;

/******************************************************************************
*                               Bremsstrahlung                                *
******************************************************************************/


// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Bremsstrahlung::Bremsstrahlung(const ParticleDef& particle_def,
                               const Medium& medium,
                               const EnergyCutSettings& cuts,
                               double multiplier,
                               bool lpm)
    : Parametrization(particle_def, medium, cuts, multiplier)
    , lorenz_(false) //TODO(mario): make it use to enable Mon 2017/09/04
    , lorenz_cut_(1e6)
    , init_lpm_effect_(true)
    , lpm_(lpm)
    , eLpm_(0)
{
}

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung& brems)
    : Parametrization(brems)
    , lorenz_(brems.lorenz_)
    , lorenz_cut_(brems.lorenz_cut_)
    , init_lpm_effect_(true)
    , lpm_(true)
    , eLpm_(brems.eLpm_)
{
}

Bremsstrahlung::~Bremsstrahlung()
{
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double Bremsstrahlung::DifferentialCrossSection(double energy, double v)
{
    double aux      =   0;
    double result   =   0;

    result = CalculateParametrization(energy, v);

    aux =   2*(components_[component_index_]->GetNucCharge())*(ME/particle_def_.mass)*RE;
    aux *=  (ALPHA/v)*aux*result;

    if(lpm_)
    {
        aux *=  lpm(energy, v);
    }

    double c2   =   pow(particle_def_.charge , 2);

    return multiplier_ * medium_->GetMolDensity()*components_[component_index_]->GetAtomInMolecule()*c2*c2*aux;
}

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Bremsstrahlung::GetIntegralLimits(double energy)
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

// ------------------------------------------------------------------------- //
double Bremsstrahlung::lpm(double energy, double v)
{
    if(init_lpm_effect_)
    {
        lpm_ = false;
        init_lpm_effect_    =   false;

        double sum      =   0;

        Integral integral_temp = Integral(IROMB,IMAXS,IPREC);

        unsigned int tmp_index = component_index_;

        for ( unsigned int i = 0; i < components_.size(); ++i)
        {
            component_index_ = i;
            Parametrization::IntegralLimits limits = GetIntegralLimits(BIGENERGY);

            sum += integral_temp.Integrate(
                limits.vMin, limits.vUp, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy, _1), 2);
            sum += integral_temp.Integrate(
                limits.vUp, limits.vMax, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy, _1), 4);
        }

        eLpm_ = ALPHA * (particle_def_.mass);
        eLpm_ *= eLpm_ / (4 * PI * ME * RE * sum);

        component_index_ = tmp_index;
        lpm_ = true;
    }

    double G, fi, xi, sp, h, s, s2, s3, ps, Gamma, Z3, Dn, s1;

    const double fi1 = 1.54954;
    const double G1  = 0.710390;
    const double G2  = 0.904912;

    Z3  =   pow((components_[component_index_]->GetNucCharge()), -1./3);

    s1  =   (components_[component_index_]->GetLogConstant())*Z3;
    Dn  =   1.54*pow((components_[component_index_]->GetAtomicNum()) , 0.27);
    s1  =   ME*Dn/((particle_def_.mass)*s1);

    s1 *= s1 * SQRT2;
    sp = sqrt(eLpm_ * v / (8 * (energy) * (1 - v)));
    h  = log(sp) / log(s1);

    if (sp < s1)
    {
        xi = 2;
    } else if (sp < 1)
    {
        xi = 1 + h - 0.08 * (1 - h) * (1 - (1 - h) * (1 - h)) / log(s1);
    } else
    {
        xi = 1;
    }

    s     = sp / sqrt(xi);
    Gamma = RE * ME / (ALPHA * (particle_def_.mass) * v);
    Gamma = 1 + 4 * PI * (medium_->GetMolDensity()) * (medium_->GetSumCharge()) * RE * pow(Gamma, 2);
    s *= Gamma;
    s2 = pow(s, 2);
    s3 = pow(s, 3);

    if (s < fi1)
    {
        fi = 1 - exp(-6 * s * (1 + (3 - PI) * s) + s3 / (0.623 + 0.796 * s + 0.658 * s2));
    } else
    {
        fi = 1 - 0.012 / pow(s2, 2);
    }

    if (s < G1)
    {
        ps = 1 - exp(-4 * s - 8 * s2 / (1 + 3.936 * s + 4.97 * s2 - 0.05 * s3 + 7.50 * pow(s2, 2)));
        G  = 3 * ps - 2 * fi;
    } else if (s < G2)
    {
        G = 36 * s2 / (36 * s2 + 1);
    } else
    {
        G = 1 - 0.022 / pow(s2, 2);
    }

    return ((xi / 3) * ((v * v) * G / (Gamma * Gamma) + 2 * (1 + (1 - v) * (1 - v)) * fi / Gamma)) /
           ((4. / 3) * (1 - v) + v * v);
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

size_t Bremsstrahlung::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    boost::hash_combine(seed, lpm_);
    boost::hash_combine(seed, lorenz_);

    return seed;
}

/******************************************************************************
*                          Specifc Parametrizations                           *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

BREMSSTRAHLUNG_IMPL(PetrukhinShestakov)
BREMSSTRAHLUNG_IMPL(KelnerKokoulinPetrukhin)
BREMSSTRAHLUNG_IMPL(CompleteScreening)
BREMSSTRAHLUNG_IMPL(AndreevBezrukovBugaev)

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

#undef BREMSSTRAHLUNG_IMPL
