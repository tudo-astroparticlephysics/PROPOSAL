
#include <cmath>
#include <stdexcept>

#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Constants.h"

#define EPAIR_PARAM_INTEGRAL_IMPL(param)                                       \
    crosssection::Epair##param::Epair##param(bool lpm)                         \
        : EpairProductionRhoIntegral(lpm)                                      \
    {                                                                          \
        hash_combine(hash, std::string(#param));                               \
    }                                                                          \
                                                                               \
    crosssection::Epair##param::Epair##param(bool lpm, const ParticleDef& p,   \
        const Medium& medium, double density_distribution)                     \
        : EpairProductionRhoIntegral(lpm, p, medium, density_distribution)     \
    {                                                                          \
        hash_combine(hash, std::string(#param));                               \
    }                                                                          \
                                                                               \
    std::unique_ptr<crosssection::Parametrization<Component>>                  \
        crosssection::Epair##param::clone() const                              \
    {                                                                          \
        using param_t                                                          \
            = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;         \
        return std::make_unique<param_t>(*this);                               \
    }

using namespace PROPOSAL;
using std::logic_error;
using std::make_tuple;

crosssection::EpairProduction::EpairProduction(bool lpm)
{
    if (lpm)
        throw std::invalid_argument("Missing particle_def and medium for Epair "
                                    "constructor with lpm=true");
}

crosssection::EpairProduction::EpairProduction(bool lpm, const ParticleDef& p,
    const Medium& medium, double density_correction)
{
    if (lpm) {
        lpm_ = std::make_shared<EpairLPM>(p, medium, density_correction);
        hash_combine(hash, lpm_->GetHash());
    }
    hash_combine(hash, p.GetHash(), medium.GetHash(), density_correction);
}

double crosssection::EpairProduction::GetLowerEnergyLim(
    const ParticleDef& p_def) const noexcept
{
    return p_def.mass + 2.f * ME;
}
crosssection::KinematicLimits crosssection::EpairProduction::GetKinematicLimits(
    const ParticleDef& p_def, const Component& comp,
    double energy) const noexcept
{
    auto lim = KinematicLimits();
    auto aux = p_def.mass / energy;
    lim.v_min = 4 * ME / energy;
    lim.v_max = 1 - 0.75 * SQRTE * aux * std::pow(comp.GetNucCharge(), 1. / 3);
    aux = 1 - 6 * aux * aux;
    lim.v_max = std::min(lim.v_max, aux);
    // limits.vMax = std::min(limits.vMax, 1 - p_def.mass / particle_energy);
    if (lim.v_max < lim.v_min)
        lim.v_max = lim.v_min;
    return lim;
}

namespace PROPOSAL {
namespace detail {
    double integrate_dedx_epair(Integral& integral,
        crosssection::EpairProduction const& param, const ParticleDef& p_def,
        const Component& comp, double energy, double v_min, double v_max)
    {
        double r1 = 0.8;
        double rUp = v_max * (1 - HALF_PRECISION);
        bool rflag = false;
        if (r1 < rUp) {
            if (2 * param.FunctionToDEdxIntegral(p_def, comp, energy, r1)
                < param.FunctionToDEdxIntegral(p_def, comp, energy, rUp)) {
                rflag = true;
            }
        }
        auto func = [&param, &p_def, &comp, energy](double v) {
            return param.FunctionToDEdxIntegral(p_def, comp, energy, v);
        };
        if (rflag) {
            if (r1 > v_max) {
                r1 = v_max;
            }
            if (r1 < v_min) {
                r1 = v_min;
            }
            auto sum = integral.Integrate(v_min, r1, func, 4);
            double r2 = std::max(1 - v_max, COMPUTER_PRECISION);
            if (r2 > 1 - r1) {
                r2 = 1 - r1;
            }
            auto func_reverse = [&, energy](double v) {
                return (1 - v)
                    * param.DifferentialCrossSection(
                        p_def, comp, energy, 1 - v);
            };
            sum += integral.Integrate(1 - v_max, r2, func_reverse, 2)
                + integral.Integrate(r2, 1 - r1, func_reverse, 4);
            return sum;
        }
        return integral.Integrate(v_min, v_max, func, 4);
    }
}
} // namespace PROPOSAL

// ------------------------------------------------------------------------- //
// Parametrization of Kelner/Kokoulin/Petrukhin
// Proc. 12th ICCR (1971), 2436
// ------------------------------------------------------------------------- //
crosssection::EpairProductionRhoIntegral::EpairProductionRhoIntegral(bool lpm)
    : crosssection::EpairProduction(lpm)
{
}

crosssection::EpairProductionRhoIntegral::EpairProductionRhoIntegral(bool lpm,
    const ParticleDef& p_def, const Medium& medium, double density_correction)
    : crosssection::EpairProduction(lpm, p_def, medium, density_correction)
{
}

// ------------------------------------------------------------------------- //
double crosssection::EpairProductionRhoIntegral::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{

    auto aux = 1 - (4 * ME) / (energy * v);
    auto aux2 = 1 - (6 * p_def.mass * p_def.mass) / (energy * energy * (1 - v));

    double rMax;
    if (aux > 0 && aux2 > 0) {
        rMax = std::sqrt(aux) * aux2;
    } else {
        rMax = 0;
    }

    aux = std::max(1 - rMax, COMPUTER_PRECISION);
    Integral integral(IROMB, IMAXS, IPREC);

    auto func = [this, &p_def, &comp, energy, v](double r) {
        return FunctionToIntegral(p_def, comp, energy, v, r);
    };

    return NA / comp.GetAtomicNum()
        * (integral.Integrate(1 - rMax, aux, func, 2)
            + integral.Integrate(aux, 1, func, 4));
}

/******************************************************************************
 *                          Specifc Parametrizations                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

EPAIR_PARAM_INTEGRAL_IMPL(KelnerKokoulinPetrukhin)
EPAIR_PARAM_INTEGRAL_IMPL(SandrockSoedingreksoRhode)
EPAIR_PARAM_INTEGRAL_IMPL(ForElectronPositron)

// ------------------------------------------------------------------------- //
double crosssection::EpairKelnerKokoulinPetrukhin::FunctionToIntegral(
    const ParticleDef& p_def, const Component& comp, double energy, double v,
    double r) const
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

    auto medium_charge = comp.GetNucCharge();
    auto medium_log_constant = comp.GetLogConstant();

    r = 1 - r; // only for integral optimization - do not forget to swap
               // integration limits!
    r2 = r * r;
    double Z3 = std::pow(medium_charge, -1. / 3);
    aux = (p_def.mass * v) / (2 * ME);
    double xi = aux * aux * (1 - r2) / (1 - v);
    double beta = (v * v) / (2 * (1 - v));

    // these are the Y_e and Y_mu expressions in the original paper
    diagram_e = (5 - r2 + 4 * beta * (1 + r2))
        / (2 * (1 + 3 * beta) * std::log(3 + 1 / xi) - r2
            - 2 * beta * (2 - r2));
    diagram_mu = (4 + r2 + 3 * beta * (1 + r2))
        / ((1 + r2) * (1.5 + 2 * beta) * std::log(3 + xi) + 1 - 1.5 * r2);

    // these arae the L_e and L_mu expressions
    aux = (1.5 * ME) / (p_def.mass * Z3);
    aux1 = (1 + xi) * (1 + diagram_e);
    aux2
        = (2 * ME * SQRTE * medium_log_constant * Z3) / (energy * v * (1 - r2));
    diagram_e = std::log((medium_log_constant * Z3 * std::sqrt(aux1))
                    / (1 + aux2 * aux1))
        - 0.5 * std::log(1 + aux * aux * aux1);
    diagram_mu = std::log((medium_log_constant / aux * Z3)
        / (1 + aux2 * (1 + xi) * (1 + diagram_mu)));

    // these are the Phi_e and Phi_mu expressions
    // if the logarithms above are below zero, the contribution of the diagram
    // is set to zero if lpm supression is taken into account, Phi_e is changed
    if (diagram_e > 0) {
        if (1 / xi < HALF_PRECISION) {
            // TODO: where does this short expression come from?
            diagram_e = (1.5 - r2 / 2 + beta * (1 + r2)) / xi * diagram_e;
        } else {
            diagram_e = (((2 + r2) * (1 + beta) + xi * (3 + r2))
                                * std::log(1 + 1 / xi)
                            + (1 - r2 - beta) / (1 + xi) - (3 + r2))
                * diagram_e;
        }
    } else {
        diagram_e = 0;
    }

    if (diagram_mu > 0) {
        diagram_mu
            = (((1 + r2) * (1 + 1.5 * beta) - (1 + 2 * beta) * (1 - r2) / xi)
                      * std::log(1 + xi)
                  + xi * (1 - r2 - beta) / (1 + xi) + (1 + 2 * beta) * (1 - r2))
            * diagram_mu;
    } else {
        diagram_mu = 0;
    }

    // Calculating the contribution of atomic electrons as scattering partner
    // Phys. Atom. Nucl. 61 (1998), 448
    if (medium_charge == 1) {
        g1 = 4.4e-5;
        g2 = 4.8e-5;
    } else {
        g1 = 1.95e-5;
        g2 = 5.3e-5;
    }

    aux = energy / p_def.mass;
    aux1 = 0.073 * std::log(aux / (1 + g1 * aux / (Z3 * Z3))) - 0.26;
    aux2 = 0.058 * std::log(aux / (1 + g2 * aux / Z3)) - 0.14;

    if (aux1 > 0 && aux2 > 0) {
        atomic_electron_contribution = aux1 / aux2;
    } else {
        atomic_electron_contribution = 0;
    }

    // combining the results
    aux = ALPHA * RE * p_def.charge;
    aux *= aux / (1.5 * PI) * 2 * medium_charge
        * (medium_charge + atomic_electron_contribution);
    aux1 = ME / p_def.mass * p_def.charge;

    aux *= (1 - v) / v * (diagram_e + aux1 * aux1 * diagram_mu);

    if (lpm_) {
        aux *= lpm_->suppression_factor(energy, v, r2, beta, xi);
    }

    if (aux < 0) {
        aux = 0;
    }

    return aux;
}

// ------------------------------------------------------------------------- //
double crosssection::EpairSandrockSoedingreksoRhode::FunctionToIntegral(
    const ParticleDef& p_def, const Component& comp, double energy, double v,
    double rho) const
{
    double m_in = p_def.mass;

    double nucl_Z = comp.GetNucCharge();
    double nucl_A = comp.GetAtomicNum();

    double rad_log = comp.GetLogConstant();

    double const_prefactor
        = 4.0 / (3.0 * PI) * nucl_Z * std::pow(ALPHA * RE * p_def.charge, 2);
    double Z13 = std::pow(nucl_Z, -1.0 / 3.0);
    double d_n = 1.54 * std::pow(nucl_A, 0.27);

    rho = 1 - rho;
    double rho2 = rho * rho;

    // --------------------------------------------------------------------- //
    // Zeta
    // --------------------------------------------------------------------- //

    double g1 = 1.95e-5;
    double g2 = 5.3e-5;

    if (nucl_Z == 1.0) {
        g1 = 4.4e-5;
        g2 = 4.8e-5;
    }

    double zeta, zeta1, zeta2;
    zeta1 = (0.073
            * std::log(energy / m_in
                / (1.0 + g1 * std::pow(nucl_Z, 2.0 / 3.0) * energy / m_in))
        - 0.26);
    zeta2 = (0.058 * std::log(energy / m_in / (1 + g2 / Z13 * energy / m_in))
        - 0.14);

    if (zeta1 > 0.0 && zeta2 > 0.0) {
        zeta = zeta1 / zeta2;
    } else {
        zeta = 0.0;
    }

    double beta = v * v / (2.0 * (1.0 - v));
    double xi = std::pow(m_in * v / (2.0 * ME), 2) * (1.0 - rho2) / (1.0 - v);

    // --------------------------------------------------------------------- //
    // Diagram e
    // --------------------------------------------------------------------- //

    double Be = ((2.0 + rho2) * (1.0 + beta) + xi * (3.0 + rho2))
            * std::log(1.0 + 1.0 / xi)
        + (1.0 - rho2 - beta) / (1.0 + xi) - (3.0 + rho2);

    double Ce2 = ((1.0 - rho2) * (1.0 + beta) + xi * (3.0 - rho2))
            * std::log(1.0 + 1.0 / xi)
        + 2.0 * (1.0 - beta - rho2) / (1.0 + xi) - (3.0 - rho2);
    double Ce1 = Be - Ce2;

    double De = ((2.0 + rho2) * (1.0 + beta) + xi * (3.0 + rho2))
            * dilog(1.0 / (1.0 + xi))
        - (2.0 + rho2) * xi * std::log(1.0 + 1.0 / xi)
        - (xi + rho2 + beta) / (1.0 + xi);

    double Le1, Le2;

    if (De / Be > 0.) {
        double Xe = std::exp(-De / Be);
        Le1 = std::log(rad_log * Z13 * std::sqrt(1.0 + xi)
                  / (Xe
                      + 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi)
                          / (energy * v * (1.0 - rho2))))
            - De / Be
            - 0.5 * std::log(Xe + std::pow(ME / m_in * d_n, 2) * (1.0 + xi));

        Le2 = std::log(rad_log * Z13 * std::exp(-1.0 / 6.0) * std::sqrt(1 + xi)
                  / (Xe
                      + 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13
                          * (1.0 + xi) / (energy * v * (1.0 - rho2))))
            - De / Be
            - 0.5
                * std::log(Xe
                    + pow(ME / m_in * d_n, 2.0) * std::exp(-1.0 / 3.0)
                        * (1.0 + xi));
    } else {
        double Xe_inv = std::exp(De / Be);
        Le1 = std::log(rad_log * Z13 * std::sqrt(1.0 + xi)
                  / (1.
                      + Xe_inv * 2.0 * ME * std::exp(0.5) * rad_log * Z13
                          * (1.0 + xi) / (energy * v * (1.0 - rho2))))
            - 0.5 * De / Be
            - 0.5
                * std::log(
                    1. + Xe_inv * std::pow(ME / m_in * d_n, 2) * (1.0 + xi));

        Le2 = std::log(rad_log * Z13 * std::exp(-1.0 / 6.0) * std::sqrt(1 + xi)
                  / (1.
                      + Xe_inv * 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13
                          * (1.0 + xi) / (energy * v * (1.0 - rho2))))
            - 0.5 * De / Be
            - 0.5
                * std::log(1.
                    + Xe_inv * std::pow(ME / m_in * d_n, 2)
                        * std::exp(-1.0 / 3.0) * (1.0 + xi));
    }

    double diagram_e = std::max(0.0,
        const_prefactor * (nucl_Z + zeta) * (1.0 - v) / v
            * (Ce1 * Le1 + Ce2 * Le2));

    // --------------------------------------------------------------------- //
    // Diagram mu
    // --------------------------------------------------------------------- //

    double Bm = ((1.0 + rho2) * (1.0 + (3.0 * beta) / 2)
                    - 1.0 / xi * (1.0 + 2.0 * beta) * (1.0 - rho2))
            * std::log(1.0 + xi)
        + xi * (1.0 - rho2 - beta) / (1.0 + xi)
        + (1.0 + 2.0 * beta) * (1.0 - rho2);

    double Cm2 = ((1.0 - beta) * (1.0 - rho2) - xi * (1.0 + rho2))
            * std::log(1.0 + xi) / xi
        - 2.0 * (1.0 - beta - rho2) / (1.0 + xi) + 1.0 - beta
        - (1.0 + beta) * rho2;
    double Cm1 = Bm - Cm2;

    double Dm = ((1.0 + rho2) * (1.0 + (3.0 * beta) / 2.0)
                    - 1.0 / xi * (1.0 + 2.0 * beta) * (1.0 - rho2))
            * dilog(xi / (1.0 + xi))
        + (1.0 + (3.0 * beta) / 2.0) * (1.0 - rho2) / xi * std::log(1.0 + xi)
        + (1.0 - rho2 - beta / 2.0 * (1.0 + rho2)
              + (1.0 - rho2) / (2.0 * xi) * beta)
            * xi / (1.0 + xi);

    double Lm1, Lm2;

    if (Dm / Bm > 0.0) {
        double Xm = std::exp(-Dm / Bm);
        Lm1 = std::log(Xm * m_in / ME * rad_log * Z13 / d_n
            / (Xm
                + 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi)
                    / (energy * v * (1.0 - rho2))));
        Lm2 = std::log(Xm * m_in / ME * rad_log * Z13 / d_n
            / (Xm
                + 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13 * (1.0 + xi)
                    / (energy * v * (1.0 - rho2))));
    } else {
        double Xm_inv = std::exp(Dm / Bm);
        Lm1 = std::log(m_in / ME * rad_log * Z13 / d_n
            / (1.0
                + 2.0 * ME * std::exp(0.5) * rad_log * Z13 * (1.0 + xi)
                    / (energy * v * (1.0 - rho2)) * Xm_inv));
        Lm2 = std::log(m_in / ME * rad_log * Z13 / d_n
            / (1.0
                + 2.0 * ME * std::exp(1.0 / 3.0) * rad_log * Z13 * (1.0 + xi)
                    / (energy * v * (1.0 - rho2)) * Xm_inv));
    }

    double diagram_mu = std::max(0.0,
        const_prefactor * (nucl_Z + zeta) * (1.0 - v) / v
            * std::pow(ME / m_in * p_def.charge , 2) * (Cm1 * Lm1 + Cm2 * Lm2));

    double aux = diagram_e + diagram_mu;

    if (lpm_) {
        aux *= lpm_->suppression_factor(energy, v, rho2, beta, xi);
    }

    if (aux < 0.0) {
        aux = 0.0;
    }

    return aux;
}

// ------------------------------------------------------------------------- //
double crosssection::EpairForElectronPositron::FunctionToIntegral(
    const ParticleDef& p_def, const Component& comp, double energy, double v,
    double rho) const
{
    // Adaptation of the direct muon pair production cross section of muons.
    // Nuclear formfactor effects have been removed as they are negligible for
    // electrons. Based on the parametrization of Kelner/Kokoulin/Petrukhin
    // Physics of Atomic Nuclei, Vol. 63, No. 9, 2000, pp. 1603-1611. Translated
    // from Yadernaya Fizika, Vol. 63, 2000, pp. 1690-1698 Original Russian Text
    // Copyright 2000 by Kel'ner, Kokoulin, Petukhin DOI: 10.1134/1.1312894

    double aux, aux1, aux2, r2, rMax, Z3, xi, beta;
    double phi, U, U_max, X, Y;
    auto medium_charge = comp.GetNucCharge();
    auto medium_log_constant = comp.GetLogConstant();
    // double medium_log_constant = 183; // According to the paper, B is set to
    // 183

    r2 = rho * rho;
    rMax = 1 - 2 * ME / (v * energy);
    Z3 = std::pow(medium_charge, -1. / 3);
    aux = v / 2;
    xi = aux * aux * (1 - r2) / (1 - v);
    beta = (v * v) / (2 * (1 - v));

    // Phi Calculation (18)
    aux = (2 + r2) * (1 + beta) + xi * (3 + r2);
    aux *= std::log(1 + 1. / xi);

    aux1 = (1 + r2) * (1 + 1.5 * beta) - 1. / xi * (1 + 2 * beta) * (1 - r2);
    aux1 *= std::log(1 + xi);
    aux2 = -1 - 3 * r2 + beta * (1 - 2 * r2);

    phi = aux + aux1 + aux2;

    // X Calculation (22)
    Y = 12 * std::sqrt(ME / energy); //(21)
    aux = medium_log_constant * Z3;
    aux1 = 2 * SQRTE * std::pow(ME, 2) * medium_log_constant * Z3 * (1 + xi)
        * (1 + Y);
    aux2 = ME * energy * v * (1 - r2);

    U = aux / (1 + aux1 / aux2);

    xi = v * v * (1 - rMax * rMax) / (4 * (1 - v));
    aux1 = 2 * SQRTE * std::pow(ME, 2) * medium_log_constant * Z3 * (1 + xi)
        * (1 + Y);
    aux2 = ME * energy * v * (1 - rMax * rMax);
    U_max = aux / (1 + aux1 / aux2);

    X = 1 + U - U_max;

    // Combine results
    aux = ALPHA * RE * p_def.charge * medium_charge;
    aux *= 2 * aux * phi * (1 - v)
        / (1.5 * PI * v); // Factor 2: Similar to factor 2 from EPairProduction,
                          // probably from symmetry in Rho

    if (X > 0) {
        aux *= std::log(X);
    } else {
        aux = 0;
    }

    if (aux < 0) {
        aux = 0;
    }

    if (lpm_) {
        aux *= lpm_->suppression_factor(energy, v, r2, beta, xi);
    }

    return aux;
}

#undef EPAIR_PARAM_INTEGRAL_IMPL

crosssection::EpairLPM::EpairLPM(
    const ParticleDef& p_def, const Medium& medium, double density_correction)
    : mass_(p_def.mass)
    , charge_(p_def.charge)
    , mol_density_(medium.GetMolDensity())
    , density_correction_(density_correction)
    , hash(0)
{
    double sum = 0.;
    auto components = medium.GetComponents();

    for (auto comp : components) {
        sum += comp.GetNucCharge() * comp.GetNucCharge()
            * std::log(3.25 * comp.GetLogConstant()
                * std::pow(comp.GetNucCharge(), -1. / 3));
    }
    eLpm_ = mass_ / (ME * RE);
    eLpm_ *= (eLpm_ * eLpm_) * ALPHA * mass_
        / (2 * PI * mol_density_ * density_correction * charge_ * charge_
            * sum);
    hash_combine(hash, mass_, std::abs(charge_), mol_density_, density_correction_);
}

// ------------------------------------------------------------------------- //
// LPM Supression by Polityko, Takahashi, Kato, Yamada, Misaki
// J. Phys. G: Nucl Part. Phys. 28 (2002) 427
// ------------------------------------------------------------------------- //

double crosssection::EpairLPM::suppression_factor(
    double energy, double v, double r2, double b, double x) const
{
    // Ternovskii functions calculated in appendix (eq. A.2)
    double A, B, C, D, E, s;
    double s36, s6, d1, d2, atan_, log1, log2;

    s = 0.25 * std::sqrt(eLpm_ / (energy * v * (1 - r2))); // eq. 29
    s6 = 6 * s;
    atan_ = s6 * (x + 1);

    if (atan_ > 1 / COMPUTER_PRECISION) {
        return 1;
    }

    s36 = 36 * s * s;
    d1 = s6 / (s6 + 1);
    d2 = s36 / (s36 + 1);
    atan_ = std::atan(atan_) - PI / 2;
    log1 = std::log((s36 * (1 + x) * (1 + x) + 1) / (s36 * x * x));
    log2 = std::log((s6 * (1 + x) + 1) / (s6 * x));
    A = 0.5 * d2 * (1 + 2 * d2 * x) * log1 - d2
        + 6 * d2 * s * (1 + ((s36 - 1) / (s36 + 1)) * x) * atan_;
    B = d1 * (1 + d1 * x) * log2 - d1;
    C = -d2 * d2 * x * log1 + d2 - (d2 * d2 * (s36 - 1) / (6 * s)) * x * atan_;
    D = d1 - d1 * d1 * x * log2;
    E = -s6 * atan_;

    return ((1 + b) * (A + (1 + r2) * B) + b * (C + (1 + r2) * D)
               + (1 - r2) * E)
        / (((2 + r2) * (1 + b) + x * (3 + r2)) * std::log(1 + 1 / x)
            + (1 - r2 - b) / (1 + x) - (3 + r2));
}
