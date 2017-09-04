
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Bremsstrahlung::Bremsstrahlung(const ParticleDef& particle_def,
                               const Medium& medium,
                               const EnergyCutSettings& cuts,
                               Definition param_def)
    : Parametrization(particle_def, medium, cuts, param_def)
    , lorenz_(false) //TODO(mario): make it use to enable Mon 2017/09/04
    , lorenz_cut_(1e6)
    , eLpm_(0)
{
}

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung& brems)
    : Parametrization(brems)
    , lorenz_(brems.lorenz_)
    , lorenz_cut_(brems.lorenz_cut_)
    , eLpm_(brems.eLpm_)
{
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Bremsstrahlung::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 0.;
    limits.vMax =
        1 - (3. / 4) * SQRTE * (particle_def_.mass / energy) * pow((current_component_->GetNucCharge()), 1. / 3);

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
double Bremsstrahlung::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * DifferentialCrossSection(energy, variable);
}

//----------------------------------------------------------------------------//
double Bremsstrahlung::FunctionToDNdxIntegral(double energy, double variable)
{
    return param_def_.multiplier * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Bremsstrahlung::DifferentialCrossSection(double energy, double v)
{
    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double s1       =   0;

    Z3  =   pow((current_component_->GetNucCharge()), -1./3);

    result = CalculateParametrization(energy, v);

    aux =   2*(current_component_->GetNucCharge())*(ME/particle_def_.mass)*RE;
    aux *=  (ALPHA/v)*aux*result;

    if(param_def_.lpm_effect_enabled)
    {
        // if(parametrization_!=ParametrizationType::BremsKelnerKokoulinPetrukhin)
        // {
        //     s1  =   (current_component_->GetLogConstant())*Z3;
        //     Dn  =   1.54*pow((current_component_->GetAtomicNum()) , 0.27);
        //     s1  =   ME*Dn/((particle_def_.mass)*s1);
        // }
        aux *=  lpm(energy, v,s1);
    }

    double c2   =   pow(particle_def_.charge , 2);

    return medium_->GetMolDensity()*current_component_->GetAtomInMolecule()*c2*c2*aux;
}

// ------------------------------------------------------------------------- //
double Bremsstrahlung::lpm(double energy, double v, double s1)
{
    if(init_lpm_effect_)
    {
        param_def_.lpm_effect_enabled = false;
        init_lpm_effect_    =   false;

        double sum      =   0;

        Integral integral_temp = Integral(IROMB,IMAXS,IPREC);

        const ComponentVec& components = medium_->GetComponents();
        const Components::Component* tmp_component = current_component_;

        for (ComponentVec::const_iterator iter = components.begin(); iter != components.end(); ++iter)
        {
            current_component_ = *iter;

            Parametrization::IntegralLimits limits = GetIntegralLimits(BIGENERGY);

            sum += integral_temp.Integrate(
                limits.vMin, limits.vUp, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy, _1), 2);
            sum += integral_temp.Integrate(
                limits.vUp, limits.vMax, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy, _1), 4);
        }

        eLpm_ = ALPHA * (particle_def_.mass);
        eLpm_ *= eLpm_ / (4 * PI * ME * RE * sum);

        current_component_ = tmp_component;
        param_def_.lpm_effect_enabled = true;
    }

    double G, fi, xi, sp, h, s, s2, s3, ps, Gamma;

    const double fi1 = 1.54954;
    const double G1  = 0.710390;
    const double G2  = 0.904912;
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
