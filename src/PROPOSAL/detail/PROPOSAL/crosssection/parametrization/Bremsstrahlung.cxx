
#include <cmath>
#include <functional>
#include <stdexcept>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/parametrization/ParamTables.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"

#define BREMSSTRAHLUNG_IMPL(param)                                             \
    crosssection::Brems##param::Brems##param(bool lpm)                         \
    {                                                                          \
        if (lpm)                                                               \
            throw std::invalid_argument(                                       \
                "Missing particle_def and medium "                             \
                "for Bremsstrahlung constructor with lpm=true");               \
        lpm_ = nullptr;                                                        \
        hash_combine(hash, std::string(#param));                               \
    }                                                                          \
                                                                               \
    crosssection::Brems##param::Brems##param(bool lpm, const ParticleDef& p,   \
        const Medium& medium, double density_correction)                       \
    {                                                                          \
        lpm_ = nullptr;                                                        \
        if (lpm) {                                                             \
            lpm_ = std::make_shared<BremsLPM>(                                 \
                p, medium, *this, density_correction);                         \
            hash_combine(hash, lpm_->GetHash());                               \
        }                                                                      \
        hash_combine(hash, std::string(#param));                               \
    }                                                                          \
                                                                               \
    std::unique_ptr<crosssection::Parametrization<Component>>                  \
        crosssection::Brems##param::clone() const                              \
    {                                                                          \
        using param_t                                                          \
            = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;         \
        return std::make_unique<param_t>(*this);                               \
    }

using namespace PROPOSAL;

crosssection::Bremsstrahlung::Bremsstrahlung()
    : lorenz_(false)
    , lorenz_cut_(1e6)
    , lpm_(nullptr)
{
}

double crosssection::Bremsstrahlung::GetLowerEnergyLim(
    const ParticleDef& p_def) const noexcept
{
    return p_def.mass;
}

crosssection::KinematicLimits crosssection::Bremsstrahlung::GetKinematicLimits(
    const ParticleDef& p_def, const Component& comp, double energy) const
{
    // The limit is taken from the Petrukhin/Shestakov Parametrization
    auto kin_lim = KinematicLimits();
    kin_lim.v_min = 0.;
    kin_lim.v_max = 1
        - 0.75 * SQRTE * (p_def.mass / energy)
            * std::pow(comp.GetNucCharge(), 1. / 3);

    if (kin_lim.v_max < 0) {
        kin_lim.v_max = 0;
    } else if (lorenz_) {
        kin_lim.v_max = std::min(kin_lim.v_max, lorenz_cut_ / energy);
    }

    // TODO: 1 - a*x is always smaller than 1 - x if a > 1
    // and 0.75*\sqrt{e}*Z^{1/3} > 1
    // so the next line will never be called, or?
    // limits.vMax = std::min(limits.vMax, 1 - p_def.mass / energy);

    return kin_lim;
}

double crosssection::Bremsstrahlung::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    // $\frac{\alpha}{v} (2 Z_{nucl} z_{particle}^2 r_e
    // \frac{m_e}{m_{particle}})^2$ is the typical Bremsstrahlung prefactor used
    // in every Parametrization

    auto result = CalculateParametrization(p_def, comp, energy, v);

    auto aux = 2 * p_def.charge * p_def.charge * (ME / p_def.mass) * RE
        * comp.GetNucCharge();
    aux *= aux * (ALPHA / v) * result;

    if (lpm_) {
        aux *= lpm_->suppression_factor(energy, v, comp);
    }

    return NA / comp.GetAtomicNum() * aux;
}

BREMSSTRAHLUNG_IMPL(PetrukhinShestakov)
BREMSSTRAHLUNG_IMPL(KelnerKokoulinPetrukhin)
BREMSSTRAHLUNG_IMPL(CompleteScreening)
BREMSSTRAHLUNG_IMPL(AndreevBezrukovBugaev)
BREMSSTRAHLUNG_IMPL(SandrockSoedingreksoRhode)

// ------------------------------------------------------------------------- //
// PetrukhinShestakov parametrization
// Canad. J. Phys. 46 (1968), 377
// ------------------------------------------------------------------------- //

double crosssection::BremsPetrukhinShestakov::CalculateParametrization(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{

    auto Z3 = std::pow(comp.GetNucCharge(), -1. / 3);
    // least momentum transferred to the nucleus (eq. 2)
    auto delta = p_def.mass * p_def.mass * v / (2 * energy * (1 - v));

    // influence of atomic form factor
    // for nuclear charge smaller 10, the nucleus is reated pointlike
    // eq. 10
    auto Fd = 189. * Z3 / ME;
    Fd = (p_def.mass) * Fd / (1 + SQRTE * delta * Fd);
    // 189 is the radiation logarithm

    // for nuclear charge greater 10, a correction for the nuclear form factor
    // is taken into account (eq.11)
    if (comp.GetNucCharge() > 10) {
        Fd *= (2. / 3) * Z3;
    }

    // eq. 3
    return ((4. / 3) * (1 - v) + v * v) * std::log(Fd);
}

// ------------------------------------------------------------------------- //
// BremsKelnerKokoulinPetrukhin parametrization
// Moscow:Preprint/MEPhI 024-95 (1995)
// ------------------------------------------------------------------------- //

double crosssection::BremsKelnerKokoulinPetrukhin::CalculateParametrization(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    double formfactor_atomic_inelastic = 0.;
    double formfactor_nuclear_inelastic = 0.;
    double result = 0.;

    // least momentum transferred to the nucleus (eq. 7)
    double delta = p_def.mass * p_def.mass * v / (2 * energy * (1 - v));

    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);
    double Dn = 1.54 * std::pow(comp.GetAtomicNum(), 0.27);
    // elastic atomic form factor (eq. 14)
    double formfactor_atomic_elastic
        = std::log(1 + ME / (delta * SQRTE * comp.GetLogConstant() * Z3));
    // elastic nuclear form factor (eq. 18)
    double formfactor_nuclear_elastic
        = std::log(Dn / (1 + delta * (Dn * SQRTE - 2) / p_def.mass));

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = (energy - p_def.mass) * (energy + p_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double maxV = ME * (energy - p_def.mass)
        / (energy * (energy - particle_momentum + ME));

    if (v < maxV) {
        // inelastic atomic contribution (eq. 26)
        formfactor_atomic_inelastic
            = std::log(p_def.mass
                  / (delta * (delta * p_def.mass / (ME * ME) + SQRTE)))
            - std::log(1 + ME / (delta * SQRTE * comp.GetBPrime() * Z3 * Z3));
    }

    if (comp.GetNucCharge() != 1) {
        // inelastic nuclear contribution (eq. 28)
        formfactor_nuclear_inelastic = formfactor_nuclear_elastic;
        // the inelastic nuclear form factor describes the scattering at single
        // nucleons in a nucleus for Hydrogen this doesn't make sense or more
        // explicit: the min required energy to excite a proton is much higher
        // than to excite a nucleus with more then just one nucleon
    }

    // eq. 2
    result = ((4. / 3) * (1 - v) + v * v)
        * (std::log(p_def.mass / delta) - 0.5 // eq.3
            - formfactor_atomic_elastic - formfactor_nuclear_elastic
            + (formfactor_nuclear_inelastic + formfactor_atomic_inelastic)
                / comp.GetNucCharge());

    return result;
}

// ------------------------------------------------------------------------- //
// CompleteScreening parametrization (by Tsai)
// Rev. Mod. Phys. 46 (1974), 815
// eq. 3.83
// ------------------------------------------------------------------------- //

double crosssection::BremsCompleteScreening::CalculateParametrization(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    (void)p_def;
    (void)energy;

    double aux = 0;
    double result = 0;
    double Lr, fZ, Lp;

    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);

    aux = ALPHA * comp.GetNucCharge();
    aux *= aux;
    fZ = aux
        * (1 / (1 + aux) + 0.20206
            + aux * (-0.0369 + aux * (0.0083 - 0.002 * aux)));

    // check rounding
    switch ((int)(comp.GetNucCharge() + 0.5)) {
    case 1: {
        Lr = 5.31;
        Lp = 6.144;
        break;
    }
    case 2: {
        Lr = 4.79;
        Lp = 5.621;
        break;
    }

    case 3: {
        Lr = 4.74;
        Lp = 5.805;
        break;
    }

    case 4: {
        Lr = 4.71;
        Lp = 5.924;
        break;
    }

    default: {
        Lr = std::log(184.15 * Z3);
        Lp = std::log(1194 * Z3 * Z3);
        break;
    }
    }

    result
        = ((4. / 3 * (1 - v) + v * v) * (comp.GetNucCharge() * (Lr - fZ) + Lp)
              + 1. / 9 * (1 - v) * (comp.GetNucCharge() + 1))
        / comp.GetNucCharge();

    return result;
}

// ------------------------------------------------------------------------- //
// AndreevBezrukovBugaev parametrization
// Phys. Atom. Nucl. 57 (1994), 2066
// ------------------------------------------------------------------------- //

double crosssection::BremsAndreevBezrukovBugaev::CalculateParametrization(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    double aux = 0;
    double result = 0;

    // least momentum transferred to the nucleus (eq. 2.2)

    double nucl_Z = comp.GetNucCharge();

    double Z3 = std::pow(nucl_Z, -1. / 3);

    double a1 = 184.15 * Z3 / (SQRTE * ME);    // eq 2.18
    double a2 = 1194 * Z3 * Z3 / (SQRTE * ME); // eq.2.19

    // calculating the contribution of elastic nuclear and atomic form factors
    // eq. 2.30
    double qc = 1.9 * MMU * Z3;
    aux = 2 * p_def.mass / qc;
    double zeta = std::sqrt(1 + aux * aux);

    double delta = p_def.mass * p_def.mass * v / (2 * energy * (1 - v));
    double x1 = a1 * delta;
    double x2 = a2 * delta;

    double aux1, aux2, d1, d2, psi1, psi2;

    if (nucl_Z == 1) {
        d1 = 0;
        d2 = 0;
    } else {
        aux1 = std::log(p_def.mass / qc);
        aux2 = 0.5 * zeta * std::log((zeta + 1) / (zeta - 1));
        d1 = aux1 + aux2;
        d2 = aux1 + 0.5 * ((3 - zeta * zeta) * aux2 + aux * aux);
    }

    // eq. 2.20 and 2.21
    aux = p_def.mass * a1;
    aux1 = std::log(aux * aux / (1 + x1 * x1));
    aux = p_def.mass * a2;
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

    if (result < 0) {
        result = 0;
    }

    return result;
}

double crosssection::BremsSandrockSoedingreksoRhode::CalculateParametrization(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    static const double a[3] = { -0.00349, 148.84, -987.531 };
    static const double b[4] = { 0.1642, 132.573, -585.361, 1407.77 };
    static const double c[6]
        = { -2.8922, -19.0156, 57.698, -63.418, 14.1166, 1.84206 };
    static const double d[6]
        = { 2134.19, 581.823, -2708.85, 4767.05, 1.52918, 0.361933 };

    double Z = comp.GetNucCharge();
    double Z13 = std::pow(Z, -1. / 3);
    double rad_log = comp.GetLogConstant();
    double rad_log_inel = comp.GetBPrime();
    double Dn = 1.54 * std::pow(comp.GetAtomicNum(), 0.27);

    double mu_qc = p_def.mass / (MMU * std::exp(1.) / Dn);
    double rho = std::sqrt(1.0 + 4.0 * mu_qc * mu_qc);

    double log_rho = std::log((rho + 1.) / (rho - 1.));
    double delta1 = std::log(mu_qc) + 0.5 * rho * log_rho;
    double delta2 = std::log(mu_qc)
        + 0.25 * (3.0 * rho - rho * rho * rho) * log_rho + 2.0 * mu_qc * mu_qc;

    // least momentum transferred to the nucleus (eq. 7)
    double delta = p_def.mass * p_def.mass * v / (2.0 * energy * (1.0 - v));

    double phi1 = std::log(rad_log * Z13 * (p_def.mass / ME)
        / (1.0 + rad_log * Z13 * exp(0.5) * delta / ME));
    double phi2 = std::log(rad_log * Z13 * exp(-1 / 6.) * (p_def.mass / ME)
        / (1.0 + rad_log * Z13 * exp(1. / 3.) * delta / ME));
    
    if (Z == 1) {
        phi1 = phi1 - delta1;
        phi2 = phi2 - delta2;
    } else {
        phi1 = phi1 - delta1 * (1. - 1. / Z);
        phi2 = phi2 - delta2 * (1. - 1. / Z);
    }

    // s_atomic
    double square_momentum = (energy - p_def.mass) * (energy + p_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double maxV = ME * (energy - p_def.mass)
        / (energy * (energy - particle_momentum + ME));

    double s_atomic = 0.0;

    if (v < maxV) {
        double s_atomic_1 = std::log(
            p_def.mass / delta / (p_def.mass * delta / (ME * ME) + SQRTE));
        double s_atomic_2
            = std::log(1. + ME / (delta * rad_log_inel * Z13 * Z13 * SQRTE));
        s_atomic = (4. / 3. * (1. - v) + v * v) * (s_atomic_1 - s_atomic_2);
    }

    // s_rad
    double s_rad;

    if (v < .0 || v > 1.0) {
        s_rad = 0.;
    } else if (v < 0.02) {
        s_rad = a[0] + a[1] * v + a[2] * v * v;
    } else if (v >= 0.02 && v < 0.1) {
        s_rad = b[0] + b[1] * v + b[2] * v * v + b[3] * v * v * v;
    } else if (v >= 0.01 && v < 0.9) {
        s_rad = c[0] + c[1] * v + c[2] * v * v;

        double tmp = std::log(1. - v);
        s_rad += c[3] * v * std::log(v) + c[4] * tmp + c[5] * tmp * tmp;
    } else {
        s_rad = d[0] + d[1] * v + d[2] * v * v;

        double tmp = std::log(1. - v);
        s_rad += d[3] * v * std::log(v) + d[4] * tmp + d[5] * tmp * tmp;
    }

    return std::max(
        ((2.0 - 2.0 * v + v * v) * phi1 - 2.0 / 3.0 * (1. - v) * phi2)
            + 1. / Z * s_atomic + 0.25 * ALPHA * phi1 * s_rad,
        0.);
}

// ------------------------------------------------------------------------- //
// EGS4 parametrization for electrons and positrons
// CompleteScreening for above 50 MeV, emperical corrections below 50 MeV
// ------------------------------------------------------------------------- //

crosssection::BremsElectronScreening::BremsElectronScreening(bool lpm)
    : interpolant_(new Interpolant(
        A_logZ, A_energies, A_correction, 2, false, false, 2, false, false))
{
    if (lpm)
        throw std::invalid_argument("Missing particle_def and medium for "
                                    "Bremsstrahlung constructor with lpm=true");
    lpm_ = nullptr;
    hash_combine(hash, std::string("electron_screening"));
}

crosssection::BremsElectronScreening::BremsElectronScreening(bool lpm,
    const ParticleDef& p_def, const Medium& medium, double density_correction)
    : interpolant_(new Interpolant(
        A_logZ, A_energies, A_correction, 2, false, false, 2, false, false))
{
    if (lpm) {
        lpm_ = std::make_shared<BremsLPM>(
            p_def, medium, *this, density_correction);
        hash_combine(hash, lpm_->GetHash());
    } else {
        lpm_ = nullptr;
    }
    hash_combine(hash, std::string("electron_screening"));
}

std::unique_ptr<crosssection::Parametrization<Component>>
crosssection::BremsElectronScreening::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

double crosssection::BremsElectronScreening::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    auto result = CalculateParametrization(p_def, comp, energy, v);

    auto aux = p_def.charge * p_def.charge * (ME / p_def.mass) * RE;
    aux *= aux * (ALPHA / v) * result;

    if (lpm_) {
        aux *= lpm_->suppression_factor(energy, v, comp);
    }

    return NA / comp.GetAtomicNum() * aux;
}

double crosssection::BremsElectronScreening::CalculateParametrization(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{

    double aux = 0;

    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);
    double logZ = std::log(comp.GetNucCharge());

    double delta = p_def.mass * p_def.mass * v / (2 * energy * (1 - v));
    double x = 136 * Z3 * 2 * delta / ME;

    // structure functions
    double phi1, phi2;
    if (x <= 1.0) {
        aux = x * x;
        phi1 = 20.867 - 3.242 * x + 0.625 * aux;
        phi2 = 20.029 - 1.930 * x - 0.086 * aux;
    } else {
        phi1 = 21.12 - 4.184 * std::log(x + 0.952);
        phi2 = phi1;
    }

    // Coulomb correction function and empirical correction factor
    double A_fac = interpolant_->InterpolateArray(logZ, energy);

    aux = ALPHA * comp.GetNucCharge();
    aux *= aux;
    auto f_c = aux
        * (1 / (1 + aux) + 0.20206
            + aux * (-0.0369 + aux * (0.0083 - 0.002 * aux)));

    // bremsstrahlung contribution with atomic electrons as target particles
    double Lr, Lp, xi;

    switch ((int)(comp.GetNucCharge() + 0.5)) {
    case 1: {
        Lr = 5.31;
        Lp = 6.144;
        break;
    }
    case 2: {
        Lr = 4.79;
        Lp = 5.621;
        break;
    }

    case 3: {
        Lr = 4.74;
        Lp = 5.805;
        break;
    }

    case 4: {
        Lr = 4.71;
        Lp = 5.924;
        break;
    }

    default: {
        Lr = std::log(184.15 * Z3);
        Lp = std::log(1194 * Z3 * Z3);
        break;
    }
    }

    xi = Lp / (Lr - f_c);

    if (energy < 50.) {
        f_c = 0;
    }

    auto result = A_fac * comp.GetNucCharge() * (comp.GetNucCharge() + xi);

    aux = (2. - 2. * v + v * v) * (phi1 - 4. / 3 * logZ - 4 * f_c)
        - 2. / 3 * (1. - v) * (phi2 - 4. / 3 * logZ - 4 * f_c);
    result *= aux;

    return result;
}

#undef BREMSSTRAHLUNG_IMPL

crosssection::BremsLPM::BremsLPM(const ParticleDef& p_def, const Medium& medium,
    const Bremsstrahlung& param, double density_correction)
    : hash(0)
    , mass_(p_def.mass)
    , mol_density_(medium.GetMolDensity())
    , mass_density_(medium.GetMassDensity())
    , sum_charge_(medium.GetSumCharge())
    , density_correction_(density_correction)
{
    hash_combine(hash, mass_, mol_density_, mass_density_, sum_charge_,
        density_correction_, param.GetHash());
    double upper_energy = 1e14;
    Integral integral_temp = Integral(IROMB, IMAXS, IPREC);
    auto components = medium.GetComponents();

    double sum = 0.;
    for (auto comp : components) {
        auto limits = param.GetKinematicLimits(p_def, comp, upper_energy);
        double contribution = 0;
        contribution += integral_temp.Integrate(
            limits.v_min, limits.v_max,
            [&param, &p_def, &comp, &upper_energy](double v) {
                return param.FunctionToDEdxIntegral(
                    p_def, comp, upper_energy, v);
            },
            2.);
        double weight_for_loss_in_medium = medium.GetSumNucleons()
            / (comp.GetAtomInMolecule() * comp.GetAtomicNum());
        sum += contribution / weight_for_loss_in_medium;
    }
    sum = sum * (mass_density_ * density_correction_);
    eLpm_ = ALPHA * mass_;
    eLpm_ *= eLpm_ / (4. * PI * ME * RE * sum);
}

double crosssection::BremsLPM::suppression_factor(
    double energy, double v, const Component& comp) const
{
    double G, fi, xi, ps, Gamma, s1;

    const double fi1 = 1.54954;
    const double G1 = 0.710390;
    const double G2 = 0.904912;
    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);

    double Dn = 1.54 * std::pow(comp.GetAtomicNum(), 0.27);
    s1 = ME * Dn / (mass_ * Z3 * comp.GetLogConstant());
    s1 = s1 * s1 * SQRT2; // TODO: is SQRT2 correct here, this factor is not in the paper?

    // Calc xi(s') from Stanev, Vankow, Streitmatter, Ellsworth, Bowen
    // Phys. Rev. D 25 (1982), 1291
    double sp = 0.125*std::sqrt(eLpm_ * v / (energy * (1 - v)));
    double h = std::log(sp) / std::log(s1);

    if (sp < s1) {
        xi = 2;
    } else if (sp < 1) {
        xi = 1 + h - 0.08 * (1 - h) * (1 - (1 - h) * (1 - h)) / std::log(s1);
    } else {
        xi = 1;
    }

    Gamma = RE * ME / (ALPHA * mass_ * v);
    Gamma = 1
        + 4 * PI * sum_charge_ * RE * Gamma * Gamma * mol_density_
            * density_correction_;
    double s = sp / std::sqrt(xi) * Gamma;
    double s2 = s * s;

    if (s < fi1) {
        // Stanev et al.,  Phys. Rev. D 25 (1982), 1291 (eq. 14d)
        fi = 1
            - std::exp(-6 * s * (1 + (3 - PI) * s)
                + s2 * s / (0.623 + 0.796 * s + 0.658 * s2));
    } else {
        fi = 1
            - 0.012 / (s2 * s2); // Migdal, Phys. Rev. 103 (1956), 1811 (eq. 48)
    }

    if (s < G1) {
        //  Stanev et al.,  Phys. Rev. D 25 (1982), 1291 (eq. 15d)
        ps = 1
            - std::exp(-4 * s
                - 8 * s2
                    / (1 + 3.936 * s + 4.97 * s2 - 0.05 * s2 * s
                        + 7.50 * s2 * s2));
        // Klein, Rev. Mod. Phys. 71 (1999), 1501 (eq. 77)
        G = 3 * ps - 2 * fi;
    } else if (s < G2) {
        G = 36 * s2 / (36 * s2 + 1);
    } else {
        G = 1
            - 0.022 / (s2 * s2); // Migdal, Phys. Rev. 103 (1956), 1811 (eq. 48)
    }

    return ((xi / 3)
               * ((v * v) * G / (Gamma * Gamma)
                   + 2 * (1 + (1 - v) * (1 - v)) * fi / Gamma))
        / ((4. / 3) * (1 - v) + v * v);
}
