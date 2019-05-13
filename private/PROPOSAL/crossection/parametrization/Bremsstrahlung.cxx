
#include <functional>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

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
    , lorenz_(false) // TODO(mario): make it use to enable Mon 2017/09/04
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
    , init_lpm_effect_(brems.init_lpm_effect_)
    , lpm_(brems.lpm_)
    , eLpm_(brems.eLpm_)
{
}

Bremsstrahlung::~Bremsstrahlung() {}

bool Bremsstrahlung::compare(const Parametrization& parametrization) const
{
    const Bremsstrahlung* bremsstrahlung = static_cast<const Bremsstrahlung*>(&parametrization);

    if (lorenz_ != bremsstrahlung->lorenz_)
        return false;
    else if (lorenz_cut_ != bremsstrahlung->lorenz_cut_)
        return false;
    else if (init_lpm_effect_ != bremsstrahlung->init_lpm_effect_)
        return false;
    else if (lpm_ != bremsstrahlung->lpm_)
        return false;
    else if (eLpm_ != bremsstrahlung->eLpm_)
        return false;
    else
        return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double Bremsstrahlung::DifferentialCrossSection(double energy, double v)
{
    // $\frac{\alpha}{v} (2 Z_{nucl} z_{particle}^2 r_e \frac{m_e}{m_{particle}})^2$
    // is the typical Bremsstrahlung prefactor used in every Parametrization

    double aux    = 0;
    double result = 0;

    result = CalculateParametrization(energy, v);

    aux = 2 * particle_def_.charge * particle_def_.charge * (ME / particle_def_.mass) * RE *
          components_[component_index_]->GetNucCharge();
    aux *= aux * (ALPHA / v) * result;

    if (lpm_)
    {
        aux *= lpm(energy, v);
    }

    return medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * aux;
}

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Bremsstrahlung::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 0.;

    // The limit is taken from the Petrukhin/Shestakov Parametrization
    limits.vMax =
        1 - 0.75 * SQRTE * (particle_def_.mass / energy) * std::pow(components_[component_index_]->GetNucCharge(), 1. / 3);

    if (limits.vMax < 0)
    {
        limits.vMax = 0;
    } else if (lorenz_)
    {
        limits.vMax = std::min(limits.vMax, lorenz_cut_ / energy);
    }

    // TODO: 1 - a*x is always smaller than 1 - x if a > 1
    // and 0.75*\sqrt{e}*Z^{1/3} > 1
    // so the next line will never be called, or?
    // limits.vMax = std::min(limits.vMax, 1 - particle_def_.mass / energy);

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    return limits;
}

// ------------------------------------------------------------------------- //
double Bremsstrahlung::lpm(double energy, double v)
{
    if (init_lpm_effect_)
    {
        lpm_             = false;
        init_lpm_effect_ = false;

        double sum = 0;

        // high energy for the calculation of the radiation length, which 
        // converges for high energie against a fixed value.
        double upper_energy = 1e14;

        Integral integral_temp = Integral(IROMB, IMAXS, IPREC);

        unsigned int tmp_index = component_index_;

        for (unsigned int i = 0; i < components_.size(); ++i)
        {
            component_index_ = i;
            Parametrization::IntegralLimits limits = GetIntegralLimits(upper_energy);

            sum += integral_temp.Integrate(
                limits.vMin, limits.vUp, std::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, upper_energy, std::placeholders::_1), 2);
            sum += integral_temp.Integrate(
                limits.vUp, limits.vMax, std::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, upper_energy, std::placeholders::_1), 4);
        }

        eLpm_ = ALPHA * (particle_def_.mass);
        eLpm_ *= eLpm_ / (4 * PI * ME * RE * sum);

        component_index_ = tmp_index;
        lpm_             = true;
    }

    double G, fi, xi, ps, Gamma, s1;

    const double fi1 = 1.54954;
    const double G1  = 0.710390;
    const double G2  = 0.904912;

    double Z3 = std::pow(components_[component_index_]->GetNucCharge(), -1. / 3);

    double Dn = 1.54 * std::pow(components_[component_index_]->GetAtomicNum(), 0.27);
    s1        = ME * Dn / (particle_def_.mass * Z3 * components_[component_index_]->GetLogConstant());
    s1 *= s1 * SQRT2;

    // Calc xi(s') from Stanev, Vankow, Streitmatter, Ellsworth, Bowen
    // Phys. Rev. D 25 (1982), 1291
    double sp = std::sqrt(eLpm_ * v / (8 * energy * (1 - v)));
    double h  = std::log(sp) / std::log(s1);

    if (sp < s1)
    {
        xi = 2;
    } else if (sp < 1)
    {
        xi = 1 + h - 0.08 * (1 - h) * (1 - (1 - h) * (1 - h)) / std::log(s1);
    } else
    {
        xi = 1;
    }

    Gamma     = RE * ME / (ALPHA * particle_def_.mass * v);
    Gamma     = 1 + 4 * PI * medium_->GetMolDensity() * medium_->GetSumCharge() * RE * Gamma * Gamma;
    double s  = sp / std::sqrt(xi) * Gamma;
    double s2 = s * s;

    if (s < fi1)
    {
        // Stanev et al.,  Phys. Rev. D 25 (1982), 1291 (eq. 14d)
        fi = 1 - std::exp(-6 * s * (1 + (3 - PI) * s) + s2 * s / (0.623 + 0.796 * s + 0.658 * s2));
    } else
    {
        fi = 1 - 0.012 / (s2 * s2); // Migdal, Phys. Rev. 103 (1956), 1811 (eq. 48)
    }

    if (s < G1)
    {
        //  Stanev et al.,  Phys. Rev. D 25 (1982), 1291 (eq. 15d)
        ps = 1 - std::exp(-4 * s - 8 * s2 / (1 + 3.936 * s + 4.97 * s2 - 0.05 * s2 * s + 7.50 * s2 * s2));
        // Klein, Rev. Mod. Phys. 71 (1999), 1501 (eq. 77)
        G = 3 * ps - 2 * fi;
    } else if (s < G2)
    {
        G = 36 * s2 / (36 * s2 + 1);
    } else
    {
        G = 1 - 0.022 / (s2 * s2); // Migdal, Phys. Rev. 103 (1956), 1811 (eq. 48)
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
    hash_combine(seed, lpm_, lorenz_);

    return seed;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void Bremsstrahlung::print(std::ostream& os) const
{
    os << "lpm enabled: " << lpm_ << '\n';
    os << "lorentz cut [MeV]: " << lorenz_cut_ << '\n';
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
BREMSSTRAHLUNG_IMPL(SandrockSoedingreksoRhode)

// ------------------------------------------------------------------------- //
// PetrukhinShestakov parametrization
// Canad. J. Phys. 46 (1968), 377
// ------------------------------------------------------------------------- //

double BremsPetrukhinShestakov::CalculateParametrization(double energy, double v)
{
    double result = 0;
    double Fd     = 0;

    double Z3 = std::pow(components_[component_index_]->GetNucCharge(), -1. / 3);
    // least momentum transferred to the nucleus (eq. 2)
    double delta = particle_def_.mass * particle_def_.mass * v / (2 * energy * (1 - v));

    // influence of atomic form factor
    // for nuclear charge smaller 10, the nucleus is reated pointlike
    // eq. 10
    Fd = 189 * Z3 / ME;
    Fd = (particle_def_.mass) * Fd / (1 + SQRTE * delta * Fd);
    // 189 is the radiation logarithm

    // for nuclear charge greater 10, a correction for the nuclear form factor
    // is taken into account (eq.11)
    if (components_[component_index_]->GetNucCharge() > 10)
    {
        Fd *= (2. / 3) * Z3;
    }

    // eq. 3
    result = ((4. / 3) * (1 - v) + v * v) * std::log(Fd);

    return result;
}

// ------------------------------------------------------------------------- //
// BremsKelnerKokoulinPetrukhin parametrization
// Moscow:Preprint/MEPhI 024-95 (1995)
// ------------------------------------------------------------------------- //

double BremsKelnerKokoulinPetrukhin::CalculateParametrization(double energy, double v)
{
    double formfactor_atomic_inelastic  = 0.;
    double formfactor_nuclear_inelastic = 0.;
    double result                       = 0.;

    // least momentum transferred to the nucleus (eq. 7)
    double delta = particle_def_.mass * particle_def_.mass * v / (2 * energy * (1 - v));

    double Z3 = std::pow(components_[component_index_]->GetNucCharge(), -1. / 3);
    double Dn = 1.54 * std::pow(components_[component_index_]->GetAtomicNum(), 0.27);
    // elastic atomic form factor (eq. 14)
    double formfactor_atomic_elastic =
        std::log(1 + ME / (delta * SQRTE * components_[component_index_]->GetLogConstant() * Z3));
    // elastic nuclear form factor (eq. 18)
    double formfactor_nuclear_elastic = std::log(Dn / (1 + delta * (Dn * SQRTE - 2) / particle_def_.mass));

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum   = (energy - particle_def_.mass) * (energy + particle_def_.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double maxV              = ME * (energy - particle_def_.mass) / (energy * (energy - particle_momentum + ME));

    if (v < maxV)
    {
        // inelastic atomic contribution (eq. 26)
        formfactor_atomic_inelastic =
            std::log(particle_def_.mass / (delta * (delta * particle_def_.mass / (ME * ME) + SQRTE))) -
            std::log(1 + ME / (delta * SQRTE * components_[component_index_]->GetBPrime() * Z3 * Z3));
    }

    if (components_[component_index_]->GetNucCharge() != 1)
    {
        // inelastic nuclear contribution (eq. 28)
        formfactor_nuclear_inelastic = formfactor_nuclear_elastic;
        // the inelastic nuclear form factor describes the scattering at single nucleons in a nucleus
        // for Hydrogen this doesn't make sense
        // or more explicit: the min required energy to excite a proton is much higher
        // than to excite a nucleus with more then just one nucleon
    }

    // eq. 2
    result = ((4. / 3) * (1 - v) + v * v) * (std::log(particle_def_.mass / delta) - 0.5 // eq.3
                                             - formfactor_atomic_elastic - formfactor_nuclear_elastic +
                                             (formfactor_nuclear_inelastic + formfactor_atomic_inelastic) /
                                                 components_[component_index_]->GetNucCharge());

    return result;
}

// ------------------------------------------------------------------------- //
// CompleteScreening parametrization (by Tsai)
// Rev. Mod. Phys. 46 (1974), 815
// eq. 3.83
// ------------------------------------------------------------------------- //

double BremsCompleteScreening::CalculateParametrization(double energy, double v)
{
    (void)energy;

    double aux    = 0;
    double result = 0;
    double Lr, fZ, Lp;

    double Z3 = std::pow(components_[component_index_]->GetNucCharge(), -1. / 3);

    aux = ALPHA * components_[component_index_]->GetNucCharge();
    aux *= aux;
    fZ = aux * (1 / (1 + aux) + 0.20206 + aux * (-0.0369 + aux * (0.0083 - 0.002 * aux)));

    // check rounding
    switch ((int)(components_[component_index_]->GetNucCharge() + 0.5))
    {
        case 1:
        {
            Lr = 5.31;
            Lp = 6.144;
            break;
        }
        case 2:
        {
            Lr = 4.79;
            Lp = 5.621;
            break;
        }

        case 3:
        {
            Lr = 4.74;
            Lp = 5.805;
            break;
        }

        case 4:
        {
            Lr = 4.71;
            Lp = 5.924;
            break;
        }

        default:
        {
            Lr = std::log(184.15 * Z3);
            Lp = std::log(1194 * Z3 * Z3);
            break;
        }
    }

    result = ((4. / 3 * (1 - v) + v * v) * (components_[component_index_]->GetNucCharge() * (Lr - fZ) + Lp) +
              1. / 9 * (1 - v) * (components_[component_index_]->GetNucCharge() + 1)) /
             components_[component_index_]->GetNucCharge();

    return result;
}

// ------------------------------------------------------------------------- //
// AndreevBezrukovBugaev parametrization
// Phys. Atom. Nucl. 57 (1994), 2066
// ------------------------------------------------------------------------- //

double BremsAndreevBezrukovBugaev::CalculateParametrization(double energy, double v)
{
    double aux    = 0;
    double result = 0;

    // least momentum transferred to the nucleus (eq. 2.2)

    double nucl_Z = components_[component_index_]->GetNucCharge();

    double Z3 = std::pow(nucl_Z, -1. / 3);

    double a1 = 184.15 * Z3 / (SQRTE * ME);    // eq 2.18
    double a2 = 1194 * Z3 * Z3 / (SQRTE * ME); // eq.2.19

    // calculating the contribution of elastic nuclear and atomic form factors
    // eq. 2.30
    double qc   = 1.9 * MMU * Z3;
    aux         = 2 * particle_def_.mass / qc;
    double zeta = std::sqrt(1 + aux * aux);

    double delta = particle_def_.mass * particle_def_.mass * v / (2 * energy * (1 - v));
    double x1    = a1 * delta;
    double x2    = a2 * delta;

    double aux1, aux2, d1, d2, psi1, psi2;

    if (nucl_Z == 1)
    {
        d1 = 0;
        d2 = 0;
    } else
    {
        aux1 = std::log(particle_def_.mass / qc);
        aux2 = 0.5 * zeta * std::log((zeta + 1) / (zeta - 1));
        d1   = aux1 + aux2;
        d2   = aux1 + 0.5 * ((3 - zeta * zeta) * aux2 + aux * aux);
    }

    // eq. 2.20 and 2.21
    aux  = particle_def_.mass * a1;
    aux1 = std::log(aux * aux / (1 + x1 * x1));
    aux  = particle_def_.mass * a2;
    aux2 = std::log(aux * aux / (1 + x2 * x2));
    psi1 = 0.5 * ((1 + aux1) + (1 + aux2) / nucl_Z);
    psi2 = 0.5 * ((2. / 3 + aux1) + (2. / 3 + aux2) / nucl_Z);

    aux1 = x1 * std::atan(1 / x1);
    aux2 = x2 * std::atan(1 / x2);
    psi1 -= aux1 + aux2 / nucl_Z;
    aux = x1 * x1;
    psi2 += 2 * aux * (1 - aux1 + 0.75 * std::log(aux / (1 + aux)));
    aux = x2 * x2;
    psi2 += 2 * aux * (1 - aux2 + 0.75 * std::log(aux / (1 + aux))) / nucl_Z;

    psi1 -= d1;
    psi2 -= d2;
    result = (2 - 2 * v + v * v) * psi1 - (2. / 3) * (1 - v) * psi2;

    if (result < 0)
    {
        result = 0;
    }

    return result;
}

double BremsSandrockSoedingreksoRhode::CalculateParametrization(double energy, double v)
{
    static const double a[3] = {-0.00349, 148.84, -987.531};
    static const double b[4] = {0.1642, 132.573, -585.361, 1407.77};
    static const double c[6] = {-2.8922, -19.0156, 57.698, -63.418, 14.1166, 1.84206};
    static const double d[6] = {2134.19, 581.823, -2708.85, 4767.05, 1.52918, 0.361933};


    double Z = components_[component_index_]->GetNucCharge();
    double Z13 = std::pow(Z, -1. / 3);
    double rad_log = components_[component_index_]->GetLogConstant();
    double rad_log_inel = components_[component_index_]->GetBPrime();
    double Dn = 1.54 * std::pow(components_[component_index_]->GetAtomicNum(), 0.27);

    double mu_qc = particle_def_.mass / (MMU * std::exp(1.) / Dn);
    double rho = std::sqrt(1.0 + 4.0*mu_qc*mu_qc);

    double log_rho = std::log((rho + 1.) / (rho - 1.));
    double delta1 = std::log(mu_qc) + 0.5 * rho * log_rho;
    double delta2 = std::log(mu_qc) + 0.25 * (3.0*rho - rho*rho*rho) * log_rho + 2.0*mu_qc*mu_qc;

    // least momentum transferred to the nucleus (eq. 7)
    double delta = particle_def_.mass * particle_def_.mass * v / (2.0 * energy * (1.0 - v));

    double phi1 = std::log(rad_log * Z13 * (particle_def_.mass / ME) / (1.0 + rad_log * Z13 * exp(0.5) * delta / ME));
    double phi2 = std::log(rad_log * Z13 * exp(-1/6.) * (particle_def_.mass / ME) / (1.0 + rad_log * Z13 * exp(1./3.) * delta / ME));
    phi1 -= delta1 * (1. - 1./Z);
    phi2 -= delta2 * (1. - 1./Z);

    // s_atomic
    double square_momentum   = (energy - particle_def_.mass) * (energy + particle_def_.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double maxV = ME * (energy - particle_def_.mass) / (energy * (energy - particle_momentum + ME));

    double s_atomic = 0.0;

    if (v < maxV)
    {
        double s_atomic_1 = std::log(particle_def_.mass / delta / ( particle_def_.mass * delta / (ME * ME) + SQRTE));
        double s_atomic_2 = std::log(1. + ME / (delta * rad_log_inel * Z13 * Z13 * SQRTE));
        s_atomic = (4./3. * (1. - v) + v*v) * (s_atomic_1 - s_atomic_2);
    }

    // s_rad
    double s_rad;

    if (v < .0 || v > 1.0)
    {
        s_rad = 0.;
    }
    else if (v < 0.02)
    {
        s_rad = a[0] + a[1] * v + a[2] * v * v;
    }
    else if (v >= 0.02 && v < 0.1)
    {
        s_rad = b[0] + b[1] * v + b[2] * v * v + b[3] * v * v * v;
    }
    else if (v >= 0.01 && v < 0.9)
    {
        s_rad = c[0] + c[1] * v + c[2] * v * v;

        double tmp = std::log(1. - v);
        s_rad += c[3] * v * std::log(v) + c[4] * tmp + c[5] * tmp * tmp;
    }
    else
    {
        s_rad = d[0] + d[1] * v + d[2] * v * v;

        double tmp = std::log(1. - v);
        s_rad += d[3] * v * std::log(v) + d[4] * tmp + d[5] * tmp * tmp;
    }

    return std::max(((2.0 - 2.0 * v + v * v) * phi1 - 2.0/3.0 * (1. - v) * phi2) + 1./Z * s_atomic + 0.25*ALPHA * phi1 * s_rad, 0.);
}

#undef BREMSSTRAHLUNG_IMPL
