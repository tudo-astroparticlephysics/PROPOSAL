
#include <cmath>

#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"

#define EPAIR_PARAM_INTEGRAL_IMPL(param)                                                                               \
    Epair##param::Epair##param(const ParticleDef& particle_def,                                                        \
                               const Medium& medium,                                                                   \
                               const EnergyCutSettings& cuts,                                                          \
                               double multiplier,                                                                      \
                               bool lpm)                                                                               \
        : EpairProductionRhoIntegral(particle_def, medium, cuts, multiplier, lpm)                                      \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Epair##param::Epair##param(const Epair##param& photo)                                                              \
        : EpairProductionRhoIntegral(photo)                                                                            \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Epair##param::~Epair##param() {}                                                                                   \
                                                                                                                       \
    const std::string Epair##param::name_ = "Epair" #param;

using namespace PROPOSAL;

/******************************************************************************
 *                              EpairProduction                               *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

EpairProduction::EpairProduction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 bool lpm)
    : Parametrization(particle_def, medium, cuts, multiplier)
    , init_lpm_effect_(true)
    , lpm_(lpm)
    , eLpm_(0)
{
}

EpairProduction::EpairProduction(const EpairProduction& epair)
    : Parametrization(epair)
    , init_lpm_effect_(epair.init_lpm_effect_)
    , lpm_(epair.lpm_)
    , eLpm_(epair.eLpm_)
{
}

EpairProduction::~EpairProduction() {}

bool EpairProduction::compare(const Parametrization& parametrization) const
{
    const EpairProduction* pairproduction = static_cast<const EpairProduction*>(&parametrization);

    if (init_lpm_effect_ != pairproduction->init_lpm_effect_)
        return false;
    else if (lpm_ != pairproduction->lpm_)
        return false;
    else if (eLpm_ != pairproduction->eLpm_)
        return false;
    else
        return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits EpairProduction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux = particle_def_.mass / energy;

    limits.vMin = 4 * ME / energy;
    limits.vMax = 1 - 0.75 * SQRTE * aux * std::pow(components_[component_index_]->GetNucCharge(), 1. / 3);

    aux         = 1 - 6 * aux * aux;
    limits.vMax = std::min(limits.vMax, aux);
    // limits.vMax = std::min(limits.vMax, 1 - particle_mass / particle_energy);

    if (limits.vMax < limits.vMin)
    {
        limits.vMax = limits.vMin;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}

// ------------------------------------------------------------------------- //
// LPM Supression by Polityko, Takahashi, Kato, Yamada, Misaki
// J. Phys. G: Nucl Part. Phys. 28 (2002) 427
// ------------------------------------------------------------------------- //

double EpairProduction::lpm(double energy, double v, double r2, double b, double x)
{
    if (init_lpm_effect_)
    {
        init_lpm_effect_ = false;
        double sum       = 0;

        for (auto component: medium_->GetComponents())
        {
            sum += component->GetNucCharge() * component->GetNucCharge() *
                   std::log(3.25 * component->GetLogConstant() * std::pow(component->GetNucCharge(), -1. / 3));
        }

        // eq. 29
        eLpm_ = particle_def_.mass / (ME * RE);
        eLpm_ *= (eLpm_ * eLpm_) * ALPHA * particle_def_.mass /
                 (2 * PI * medium_->GetMolDensity() * particle_def_.charge * particle_def_.charge * sum);
    }

    // Ternovskii functions calculated in appendix (eq. A.2)
    double A, B, C, D, E, s;
    double s36, s6, d1, d2, atan_, log1, log2;

    s     = 0.25 * std::sqrt(eLpm_ / (energy * v * (1 - r2))); // eq. 29
    s6    = 6 * s;
    atan_ = s6 * (x + 1);

    if (atan_ > 1 / COMPUTER_PRECISION)
    {
        return 1;
    }

    s36   = 36 * s * s;
    d1    = s6 / (s6 + 1);
    d2    = s36 / (s36 + 1);
    atan_ = std::atan(atan_) - PI / 2;
    log1  = std::log((s36 * (1 + x) * (1 + x) + 1) / (s36 * x * x));
    log2  = std::log((s6 * (1 + x) + 1) / (s6 * x));
    A     = 0.5 * d2 * (1 + 2 * d2 * x) * log1 - d2 + 6 * d2 * s * (1 + ((s36 - 1) / (s36 + 1)) * x) * atan_;
    B     = d1 * (1 + d1 * x) * log2 - d1;
    C     = -d2 * d2 * x * log1 + d2 - (d2 * d2 * (s36 - 1) / (6 * s)) * x * atan_;
    D     = d1 - d1 * d1 * x * log2;
    E     = -s6 * atan_;

    return ((1 + b) * (A + (1 + r2) * B) + b * (C + (1 + r2) * D) + (1 - r2) * E) /
           (((2 + r2) * (1 + b) + x * (3 + r2)) * std::log(1 + 1 / x) + (1 - r2 - b) / (1 + x) - (3 + r2));
}

// ------------------------------------------------------------------------- //
// Parametrization of Kelner/Kokoulin/Petrukhin
// Proc. 12th ICCR (1971), 2436
// ------------------------------------------------------------------------- //

/******************************************************************************
 *                          EpairProduction Integral                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
EpairProductionRhoIntegral::EpairProductionRhoIntegral(const ParticleDef& particle_def,
                                                       const Medium& medium,
                                                       const EnergyCutSettings& cuts,
                                                       double multiplier,
                                                       bool lpm)
    : EpairProduction(particle_def, medium, cuts, multiplier, lpm)
    , integral_(IROMB, IMAXS, IPREC)
{
}

// ------------------------------------------------------------------------- //
EpairProductionRhoIntegral::EpairProductionRhoIntegral(const EpairProductionRhoIntegral& epair)
    : EpairProduction(epair)
    , integral_(epair.integral_)
{
}

// ------------------------------------------------------------------------- //
EpairProductionRhoIntegral::~EpairProductionRhoIntegral() {}

// ------------------------------------------------------------------------- //
bool EpairProductionRhoIntegral::compare(const Parametrization& parametrization) const
{
    const EpairProductionRhoIntegral* epair = static_cast<const EpairProductionRhoIntegral*>(&parametrization);

    if (integral_ != epair->integral_)
        return false;
    else
        return EpairProduction::compare(parametrization);
}

// ------------------------------------------------------------------------- //
double EpairProductionRhoIntegral::DifferentialCrossSection(double energy, double v)
{
    double rMax, aux, aux2;

    aux  = 1 - (4 * ME) / (energy * v);
    aux2 = 1 - (6 * particle_def_.mass * particle_def_.mass) / (energy * energy * (1 - v));

    if (aux > 0 && aux2 > 0)
    {
        rMax = std::sqrt(aux) * aux2;
    } else
    {
        rMax = 0;
    }

    aux = std::max(1 - rMax, COMPUTER_PRECISION);

    return medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() *
           particle_def_.charge * particle_def_.charge *
           (integral_.Integrate(
                1 - rMax, aux, std::bind(&EpairProductionRhoIntegral::FunctionToIntegral, this, energy, v, std::placeholders::_1), 2) +
            integral_.Integrate(
                aux, 1, std::bind(&EpairProductionRhoIntegral::FunctionToIntegral, this, energy, v, std::placeholders::_1), 4));
}

// ------------------------------------------------------------------------- //
void EpairProductionRhoIntegral::print(std::ostream& os) const
{
    os << "lpm enabled: " << lpm_ << '\n';
}

// ------------------------------------------------------------------------- //
size_t EpairProductionRhoIntegral::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed, lpm_);

    return seed;
}

/******************************************************************************
 *                          Specifc Parametrizations                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

EPAIR_PARAM_INTEGRAL_IMPL(KelnerKokoulinPetrukhin)
EPAIR_PARAM_INTEGRAL_IMPL(SandrockSoedingreksoRhode)

// ------------------------------------------------------------------------- //
double EpairKelnerKokoulinPetrukhin::FunctionToIntegral(double energy, double v, double r)
{
    // Parametrization of Kelner/Kokoulin/Petrukhin
    // Proc. 12th ICCR (1971), 2436
    //
    // there are two pair production diagrams taking into account here
    // where an electron (or positron) couples to the nucleus (marked with xx_e)
    // and where the muon couples to the nucleus (marked with xx_mu)
    //
    // an additional contribution comes from the scattering on atomic electrons
    // this is described below

    double g1, g2;
    double aux, aux1, aux2, r2;
    double diagram_e, diagram_mu, atomic_electron_contribution;
    double medium_charge       = components_[component_index_]->GetNucCharge();
    double medium_log_constant = components_[component_index_]->GetLogConstant();

    r           = 1 - r; // only for integral optimization - do not forget to swap integration limits!
    r2          = r * r;
    double Z3   = std::pow(medium_charge, -1. / 3);
    aux         = (particle_def_.mass * v) / (2 * ME);
    double xi   = aux * aux * (1 - r2) / (1 - v);
    double beta = (v * v) / (2 * (1 - v));

    // these are the Y_e and Y_mu expressions in the original paper
    diagram_e  = (5 - r2 + 4 * beta * (1 + r2)) / (2 * (1 + 3 * beta) * std::log(3 + 1 / xi) - r2 - 2 * beta * (2 - r2));
    diagram_mu = (4 + r2 + 3 * beta * (1 + r2)) / ((1 + r2) * (1.5 + 2 * beta) * std::log(3 + xi) + 1 - 1.5 * r2);

    // these arae the L_e and L_mu expressions
    aux        = (1.5 * ME) / (particle_def_.mass * Z3);
    aux1       = (1 + xi) * (1 + diagram_e);
    aux2       = (2 * ME * SQRTE * medium_log_constant * Z3) / (energy * v * (1 - r2));
    diagram_e  = std::log((medium_log_constant * Z3 * std::sqrt(aux1)) / (1 + aux2 * aux1)) - 0.5 * std::log(1 + aux * aux * aux1);
    diagram_mu = std::log((medium_log_constant / aux * Z3) / (1 + aux2 * (1 + xi) * (1 + diagram_mu)));

    // these are the Phi_e and Phi_mu expressions
    // if the logarithms above are below zero, the contribution of the diagram is set to zero
    // if lpm supression is taken into account, Phi_e is changed
    if (diagram_e > 0)
    {
        if (1 / xi < HALF_PRECISION)
        {
            // TODO: where does this short expression come from?
            diagram_e = (1.5 - r2 / 2 + beta * (1 + r2)) / xi * diagram_e;
        } else
        {
            diagram_e =
                (((2 + r2) * (1 + beta) + xi * (3 + r2)) * std::log(1 + 1 / xi) + (1 - r2 - beta) / (1 + xi) - (3 + r2)) *
                diagram_e;
        }
    } else
    {
        diagram_e = 0;
    }

    if (diagram_mu > 0)
    {
        diagram_mu = (((1 + r2) * (1 + 1.5 * beta) - (1 + 2 * beta) * (1 - r2) / xi) * std::log(1 + xi) +
                      xi * (1 - r2 - beta) / (1 + xi) + (1 + 2 * beta) * (1 - r2)) *
                     diagram_mu;
    } else
    {
        diagram_mu = 0;
    }

    // Calculating the contribution of atomic electrons as scattering partner
    // Phys. Atom. Nucl. 61 (1998), 448
    if (medium_charge == 1)
    {
        g1 = 4.4e-5;
        g2 = 4.8e-5;
    } else
    {
        g1 = 1.95e-5;
        g2 = 5.3e-5;
    }

    aux  = energy / particle_def_.mass;
    aux1 = 0.073 * std::log(aux / (1 + g1 * aux / (Z3 * Z3))) - 0.26;
    aux2 = 0.058 * std::log(aux / (1 + g2 * aux / Z3)) - 0.14;

    if (aux1 > 0 && aux2 > 0)
    {
        atomic_electron_contribution = aux1 / aux2;
    } else
    {
        atomic_electron_contribution = 0;
    }

    // combining the results
    aux = ALPHA * RE * particle_def_.charge;
    ;
    aux *= aux / (1.5 * PI) * 2 * medium_charge * (medium_charge + atomic_electron_contribution);
    aux1 = ME / particle_def_.mass * particle_def_.charge;
    ;
    aux *= (1 - v) / v * (diagram_e + aux1 * aux1 * diagram_mu);

    if (lpm_)
    {
        aux *= lpm(energy, v, r2, beta, xi);
    }

    if (aux < 0)
    {
        aux = 0;
    }

    return aux;
}

// ------------------------------------------------------------------------- //
double EpairSandrockSoedingreksoRhode::FunctionToIntegral(double energy, double v, double rho)
{
    double m_in = particle_def_.mass;

    double nucl_Z = components_[component_index_]->GetNucCharge();
    double nucl_A = components_[component_index_]->GetAtomicNum();

    double rad_log = components_[component_index_]->GetLogConstant();

    double const_prefactor = 4.0 / (3.0 * PI) * nucl_Z * std::pow(ALPHA * RE, 2.0);
    double Z13             = std::pow(nucl_Z, -1.0 / 3.0);
    double d_n             = 1.54 * std::pow(nucl_A, 0.27);

    rho         = 1 - rho;
    double rho2 = rho * rho;

    // --------------------------------------------------------------------- //
    // Zeta
    // --------------------------------------------------------------------- //

    double g1 = 1.95e-5;
    double g2 = 5.3e-5;

    if (nucl_Z == 1.0)
    {
        g1 = 4.4e-5;
        g2 = 4.8e-5;
    }

    double zeta, zeta1, zeta2;
    zeta1 = (0.073 * std::log(energy / m_in / (1.0 + g1 * std::pow(nucl_Z, 2.0 / 3.0) * energy / m_in)) - 0.26);
    zeta2 = (0.058 * std::log(energy / m_in / (1 + g2 / Z13 * energy / m_in)) - 0.14);

    if (zeta1 > 0.0 && zeta2 > 0.0)
    {
        zeta = zeta1 / zeta2;
    } else
    {
        zeta = 0.0;
    }

    double beta = v * v / (2.0 * (1.0 - v));
    double xi   = std::pow(m_in * v / (2.0 * ME), 2.0) * (1.0 - rho2) / (1.0 - v);

    // --------------------------------------------------------------------- //
    // Diagram e
    // --------------------------------------------------------------------- //

    double Be = ((2.0 + rho2) * (1.0 + beta) + xi * (3.0 + rho2)) * std::log(1.0 + 1.0 / xi) +
                (1.0 - rho2 - beta) / (1.0 + xi) - (3.0 + rho2);

    double Ce2 = ((1.0 - rho2) * (1.0 + beta) + xi * (3.0 - rho2)) * std::log(1.0 + 1.0 / xi) +
                 2.0 * (1.0 - beta - rho2) / (1.0 + xi) - (3.0 - rho2);
    double Ce1 = Be - Ce2;

    double De = ((2.0 + rho2) * (1.0 + beta) + xi * (3.0 + rho2)) * dilog(1.0 / (1.0 + xi)) -
                (2.0 + rho2) * xi * std::log(1.0 + 1.0 / xi) - (xi + rho2 + beta) / (1.0 + xi);

    double Le1, Le2;

    if (De / Be > 0.)
    {
        double Xe = std::exp(-De / Be);
        Le1 = std::log(rad_log * Z13 * std::sqrt(1.0 + xi) /
                         (Xe + 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2)))) -
                     De / Be - 0.5 * std::log(Xe + std::pow(ME / m_in * d_n, 2.0) * (1.0 + xi));

        Le2 = std::log(rad_log * Z13 * std::exp(-1.0 / 6.0) * std::sqrt(1 + xi) /
                         (Xe + 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2)))) -
                     De / Be - 0.5 * std::log(Xe + pow(ME / m_in * d_n, 2.0) * std::exp(-1.0 / 3.0) * (1.0 + xi));
    }
    else
    {
        double Xe_inv = std::exp(De / Be);
        Le1 = std::log(rad_log * Z13 * std::sqrt(1.0 + xi) /
                         (1. + Xe_inv * 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2)))) -
                     0.5 * De / Be - 0.5 * std::log(1. + Xe_inv * std::pow(ME / m_in * d_n, 2.0) * (1.0 + xi));

        Le2 = std::log(rad_log * Z13 * std::exp(-1.0 / 6.0) * std::sqrt(1 + xi) /
                         (1. + Xe_inv * 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2)))) -
                     0.5 * De / Be - 0.5 * std::log(1. + Xe_inv * std::pow(ME / m_in * d_n, 2.0) * std::exp(-1.0 / 3.0) * (1.0 + xi));
    }

    double diagram_e = std::max(0.0, const_prefactor * (nucl_Z + zeta) * (1.0 - v) / v * (Ce1 * Le1 + Ce2 * Le2));

    // --------------------------------------------------------------------- //
    // Diagram mu
    // --------------------------------------------------------------------- //

    double Bm =
        ((1.0 + rho2) * (1.0 + (3.0 * beta) / 2) - 1.0 / xi * (1.0 + 2.0 * beta) * (1.0 - rho2)) * std::log(1.0 + xi) +
        xi * (1.0 - rho2 - beta) / (1.0 + xi) + (1.0 + 2.0 * beta) * (1.0 - rho2);

    double Cm2 = ((1.0 - beta) * (1.0 - rho2) - xi * (1.0 + rho2)) * std::log(1.0 + xi) / xi -
                 2.0 * (1.0 - beta - rho2) / (1.0 + xi) + 1.0 - beta - (1.0 + beta) * rho2;
    double Cm1 = Bm - Cm2;

    double Dm = ((1.0 + rho2) * (1.0 + (3.0 * beta) / 2.0) - 1.0 / xi * (1.0 + 2.0 * beta) * (1.0 - rho2)) *
                    dilog(xi / (1.0 + xi)) +
                (1.0 + (3.0 * beta) / 2.0) * (1.0 - rho2) / xi * std::log(1.0 + xi) +
                (1.0 - rho2 - beta / 2.0 * (1.0 + rho2) + (1.0 - rho2) / (2.0 * xi) * beta) * xi / (1.0 + xi);

    double Lm1, Lm2;

    if (Dm / Bm > 0.0)
    {
        double Xm = std::exp(-Dm / Bm);
        Lm1 = std::log(Xm * m_in / ME * rad_log * Z13 / d_n /
                  (Xm + 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2))));
        Lm2 = std::log(Xm * m_in / ME * rad_log * Z13 / d_n /
                  (Xm + 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2))));
    } else
    {
        double Xm_inv = std::exp(Dm / Bm);
        Lm1 = std::log(m_in / ME * rad_log * Z13 / d_n /
                  (1.0 + 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2)) * Xm_inv));
        Lm2 = std::log(
            m_in / ME * rad_log * Z13 / d_n /
            (1.0 + 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13 * (1.0 + xi) / (energy * v * (1.0 - rho2)) * Xm_inv));
    }

    double diagram_mu = std::max(
        0.0, const_prefactor * (nucl_Z + zeta) * (1.0 - v) / v * std::pow(ME / m_in, 2.0) * (Cm1 * Lm1 + Cm2 * Lm2));

    double aux = diagram_e + diagram_mu;

    if (lpm_)
    {
        aux *= lpm(energy, v, rho2, beta, xi);
    }

    if (aux < 0.0)
    {
        aux = 0.0;
    }

    return aux;
}

#undef EPAIR_PARAM_INTEGRAL_IMPL
