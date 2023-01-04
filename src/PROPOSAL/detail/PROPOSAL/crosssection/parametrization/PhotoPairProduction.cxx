
#include <cmath>

#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/crosssection/parametrization/ParamTables.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/Integral.h"

using namespace PROPOSAL;
using std::make_tuple;

crosssection::PhotoPairProduction::PhotoPairProduction()
    : lpm_(nullptr)
    , density_correction_(1.0)
{}

double crosssection::PhotoPairProduction::GetLowerEnergyLim(
    const ParticleDef&) const noexcept
{
    return 2. * ME;
}

crosssection::KinematicLimits
crosssection::PhotoPairProduction::GetKinematicLimits(
    const ParticleDef&, const Component&, double energy) const noexcept
{
    // x is the integration variable here
    KinematicLimits lim;
    if (energy <= 2. * ME) {
        lim.v_min = 0.5;
        lim.v_max = 0.5;
    } else {
        lim.v_min = ME / energy;
        lim.v_max = 1. - ME / energy;
    }
    return lim;
}

crosssection::PhotoPairTsai::PhotoPairTsai(bool lpm)
{
    if (lpm)
        throw std::invalid_argument("Missing particle_def and medium for "
                                    "PhotoPairProduction constructor with "
                                    "lpm=true");
    lpm_ = nullptr;
    hash_combine(hash, std::string("tsai"));
}

crosssection::PhotoPairTsai::PhotoPairTsai(
        bool lpm, const ParticleDef& p_def, const Medium& medium,
        double density_correction)
{
    if (lpm) {
        lpm_ = std::make_shared<PhotoPairLPM>(p_def, medium, *this);
        hash_combine(hash, density_correction, lpm_->GetHash());
        density_correction_ = density_correction;
    } else {
        lpm_ = nullptr;
    }
    hash_combine(hash, std::string("tsai"));
}

std::unique_ptr<crosssection::Parametrization<Component>>
crosssection::PhotoPairTsai::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

std::unique_ptr<crosssection::Parametrization<Component>>
crosssection::PhotoPairKochMotz::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

crosssection::PhotoPairKochMotz::PhotoPairKochMotz(bool lpm)
        : interpolant_(new Interpolant(photopair_KM_Z, photopair_KM_energies,
                                     photopair_KM_cross, 2, false, false,
                                     2, false, false))
{
    if (lpm)
        throw std::invalid_argument("Missing particle_def and medium for "
                                    "PhotoPairProduction constructor with "
                                    "lpm=true");
    lpm_ = nullptr;
    hash_combine(hash, std::string("kochmotz"));
}

crosssection::PhotoPairKochMotz::PhotoPairKochMotz(
        bool lpm, const ParticleDef& p_def, const Medium& medium,
        double density_correction)
    : interpolant_(new Interpolant(photopair_KM_Z, photopair_KM_energies,
                                   photopair_KM_cross, 2, false, false,
                                   2, false, false))
{
    if (lpm) {
        lpm_ = std::make_shared<PhotoPairLPM>(p_def, medium, *this);
        hash_combine(hash, density_correction, lpm_->GetHash());
        density_correction_ = density_correction;
    } else {
        lpm_ = nullptr;
    }
    hash_combine(hash, std::string("kochmotz"));
}

double crosssection::PhotoPairKochMotz::DifferentialCrossSection(
        const ParticleDef& p, const Component& comp, double energy, double v) const
{
    if (energy > 50)
        return DifferentialCrossSectionWithoutA(p, comp, energy, v);

    // Correction factor A for low energies given by the Storm and Israel data
    // (doi.org/10.1016/S0092-640X(70)80017-1)
    auto limits = PhotoPairProduction::GetKinematicLimits(p, comp, energy);
    Integral i;
    auto integrand = [this, &p, &comp, energy](double v) {
        return this->DifferentialCrossSectionWithoutA(p, comp, energy, v) / NA * comp.GetAtomicNum() / (1e-24);
    };
    auto dNdx_nocorrection = i.Integrate(limits.v_min, limits.v_max, integrand, 3);
    auto A = interpolant_->InterpolateArray(comp.GetNucCharge(), energy) / dNdx_nocorrection;
    return A * DifferentialCrossSectionWithoutA(p, comp, energy, v);
}

double crosssection::PhotoPairKochMotz::DifferentialCrossSectionWithoutA(
        const ParticleDef&, const Component& comp, double energy, double v) const
{
    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);
    double logZ = std::log(comp.GetNucCharge());

    double delta = ME * ME / (2 * energy * v * (1 - v));
    double x = 136 * Z3 * 2 * delta / ME;

    // structure functions
    double phi1, phi2, aux;
    if (x <= 1.0) {
        aux = x * x;
        phi1 = 20.867 - 3.242 * x + 0.625 * aux;
        phi2 = 20.029 - 1.930 * x - 0.086 * aux;
    } else {
        phi1 = 21.12 - 4.184 * std::log(x + 0.952);
        phi2 = phi1;
    }

    aux = ALPHA * comp.GetNucCharge();
    aux *= aux;
    auto f_c = aux
               * (1 / (1 + aux) + 0.20206
                  + aux * (-0.0369 + aux * (0.0083 - 0.002 * aux)));

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

    if (energy < 50)
        f_c = 0;

    auto result = comp.GetNucCharge() * (comp.GetNucCharge() + xi);

    result *= RE * RE * ALPHA * ( (2 * v * v - 2 * v + 1) * (phi1 - 4./3. * logZ - 4. * f_c)
            + 2./3 * v * (1. - v) * (phi2 - 4./3. * logZ - 4 * f_c));

    if (lpm_) {
        result *= lpm_->suppression_factor(energy, v, comp, density_correction_);
    }

    return result * NA / comp.GetAtomicNum();

}

double crosssection::PhotoPairTsai::DifferentialCrossSection(
    const ParticleDef&, const Component& comp, double energy, double x) const
{
    // Pair production and bremsstrahlung of chraged leptons, Yung-Su Tsai,
    // Review of Modern Physics, Vol. 46, No. 4, October 1974
    // see formula (3.9)

    double Phi1, Phi2, Psi1, Psi2;
    double Z = comp.GetNucCharge();
    double logZ = std::log(comp.GetNucCharge());
    double eta;
    double k = energy;
    double delta = std::pow(ME, 2.) / (2. * k * x * (1. - x)); // (3.20);

    if (Z <= 2.5) {
        // Hydrogen or Helium
        if (Z <= 1.5) {
            eta = 1;
        } else {
            eta = 1.6875;
        }

        double C = delta / (2. * ALPHA * ME * eta);

        Phi1 = 4. / 3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA))
            + 13. / 3. - 2. * std::log(1. + std::pow(C, 2.))
            - (13. / 2.) * C * std::atan(1. / C)
            + (1. / 6.) / (1 + std::pow(C, -2.)); // (3.25), with erratum
        Phi2 = 4. / 3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA))
            + 11. / 3. - 2. * std::log(1 + std::pow(C, 2.))
            + 25. * std::pow(C, 2.) * (1. - C * std::atan(1. / C))
            - 14. * std::pow(C, 2.) * std::log(1 + std::pow(C, -2.)); // (3.26)
        Psi1 = 8. / 3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA))
            + 23. / 3. - 2. * std::log(1 + std::pow(C, 2.))
            - 17.5 * C * std::atan(1. / C)
            + 8. * std::pow(C, 2.) * std::log(1 + std::pow(C, -2.))
            - 1. / 6. / (1 + std::pow(C, -2.)); // (3.27)
        Psi2 = 8. / 3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA))
            + 21. / 3. - 2. * std::log(1 + std::pow(C, 2.))
            - 105. * std::pow(C, 2.) * (1. - C * std::atan(1. / C))
            + 50. * std::pow(C, 2.) * std::log(1. + std::pow(C, -2.))
            - 24. * std::pow(C, 2.)
                * (-std::log(std::pow(C, 2.)) * std::log(1 + std::pow(C, -2.))
                    + dilog(1. + std::pow(C, -2.)) - dilog(1)); // (3.28)
    } else if (Z <= 4.5) {
        // Lithium or Berylium
        double a, b, ap, bp;
        if (Z <= 3.5) {
            a = 100. * std::pow(Z, -1. / 3.) / ME;
            ap = 418.6 * std::pow(Z, -2. / 3.) / ME;
        } else {
            a = 106. * std::pow(Z, -1. / 3.) / ME;
            ap = 571.4 * std::pow(Z, -2. / 3.) / ME;
        }
        b = delta * a;
        bp = delta * ap;

        Phi1 = 2.
                * (1.
                    + std::log(std::pow(a, 2.) * std::pow(Z, 2. / 3.)
                        * std::pow(ME, 2.)))
            - 2. * std::log(1. + std::pow(b, 2.))
            - 4. * b * std::atan(1. / b); // (3.46)
        Phi2 = 2.
                * (2. / 3.
                    + std::log(std::pow(a, 2.) * std::pow(Z, 2. / 3.)
                        * std::pow(ME, 2.)))
            - 2. * std::log(1. + std::pow(b, 2.))
            + 8. * std::pow(b, 2.)
                * (1 - b * std::atan(1. / b)
                    - 0.75 * std::log(1 + std::pow(b, -2.))); // (3.47)
        Psi1 = 2.
                * (1.
                    + std::log(std::pow(ap, 2.) * std::pow(Z, 4. / 3.)
                        * std::pow(ME, 2.)))
            - 2. * std::log(1 + std::pow(bp, 2.))
            - 4. * bp * std::atan(1. / b); // (3.48)
        Psi2 = 2.
                * (2. / 3.
                    + std::log(std::pow(ap, 2.) * std::pow(Z, 4. / 3.)
                        * std::pow(ME, 2.)))
            - 2. * std::log(1. + std::pow(bp, 2.))
            + 8. * std::pow(bp, 2.)
                * (1. - b * std::atan(1. / bp)
                    - 0.75 * std::log(1 + std::pow(bp, -2.))); // (3.49)
    } else {
        // Heavier elements
        double gamma = 200. * delta / (ME * std::pow(Z, 1. / 3.));   // (3.30)
        double epsilon = 200. * delta / (ME * std::pow(Z, 2. / 3.)); // (3.31)

        Phi1 = 20.863 - 2. * std::log(1. + std::pow(0.55846 * gamma, 2.))
            - 4.
                * (1. - 0.6 * std::exp(-0.9 * gamma)
                    - 0.4 * std::exp(-1.5 * gamma)); // (3.38)
        Phi2 = Phi1
            - 2. / 3. * 1.
                / (1. + 6.5 * gamma + 6. * std::pow(gamma, 2.)); // (3.39)
        Psi1 = 28.340 - 2. * std::log(1. + std::pow(3.621 * epsilon, 2.))
            - 4.
                * (1. - 0.7 * std::exp(-8. * epsilon)
                    - 0.3 * std::exp(-29.2 * epsilon)); // (3.40)
        Psi2 = Psi1
            - 2. / 3. * 1.
                / (1. + 40. * epsilon + 400. * std::pow(epsilon, 2.)); // (3.41)
    }

    double z = std::pow(Z / 137., 2.);
    double f = 1.202 * z - 1.0369 * std::pow(z, 2.)
        + 1.008 * std::pow(z, 3.) / (1. + z); // (3.3)

    double aux = (4. / 3. * std::pow(x, 2.) - 4. / 3. * x + 1.)
            * (std::pow(Z, 2.) * (Phi1 - 4. / 3. * logZ - 4. * f)
                + Z * (Psi1 - 8. / 3. * logZ))
        - 2. / 3. * x * (1. - x)
            * (std::pow(Z, 2.) * (Phi1 - Phi2) + Z * (Psi1 - Psi2)); // (3.9)

    aux *= ALPHA * std::pow(RE, 2.) / k; // (3.9)

    double p = std::sqrt(
        std::pow(x * k, 2.) - std::pow(ME, 2.)); // electron momentum
    aux *= x * std::pow(k, 2.) / p; // conversion from differential cross
                                    // section in electron momentum to x

    if (lpm_) {
        aux *= lpm_->suppression_factor(energy, x, comp, density_correction_);
    }

    return std::max(NA / comp.GetAtomicNum() * aux,
        0.); // TODO what are the real factors here, those are just guesses
}

crosssection::PhotoPairLPM::PhotoPairLPM(const ParticleDef& p_def, const Medium& medium,
                                         const PhotoPairProduction& param)
    : hash(0)
    , mol_density_(medium.GetMolDensity())
    , mass_density_(medium.GetMassDensity())
    , sum_charge_(medium.GetSumCharge())
{
    hash_combine(hash, mol_density_, mass_density_, sum_charge_, param.GetHash());
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
                return param.DifferentialCrossSection(
                    p_def, comp, upper_energy, v);
            },
            2.);
        double weight_for_loss_in_medium = medium.GetSumNucleons()
            / (comp.GetAtomInMolecule() * comp.GetAtomicNum());
        sum += contribution / weight_for_loss_in_medium;
    }
    sum = 9./7. * sum * mass_density_;
    eLpm_ = ALPHA * ME;
    eLpm_ *= 2 * eLpm_ / (PI * ME * RE * sum);
}

double crosssection::PhotoPairLPM::suppression_factor(
        double energy, double x, const Component& comp,
        double density_correction) const
{
    // taken from crosssection::BremsLPM::suppression_factor with appropriate modifications
    double G, fi, xi, ps, Gamma, s1;

    const double fi1 = 1.54954;
    const double G1 = 0.710390;
    const double G2 = 0.904912;
    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);

    // electrons are much lighter than muons, therefore no nuclear formfactor correction
    s1 = 1 / (Z3 * comp.GetLogConstant()); // PRD 25 (1982), 1291, Eq. (17)
    s1 = s1 * s1 * SQRT2; // TODO: is SQRT2 correct here, this factor is not in the paper?

    // Calc xi(s') from Stanev, Vankow, Streitmatter, Ellsworth, Bowen
    // Phys. Rev. D 25 (1982), 1291, Eq. (19)
    double sp = 0.125 * std::sqrt(eLpm_ / (density_correction * energy * x * (1 - x)));
    double h = std::log(sp) / std::log(s1);

    if (sp < s1) {
        xi = 2;
    } else if (sp < 1) {
        xi = 1 + h - 0.08 * (1 - h) * (1 - (1 - h) * (1 - h)) / std::log(s1);
    } else {
        xi = 1;
    }

    Gamma = RE * ME / (ALPHA * ME * x);
    Gamma = 1
        + 4 * PI * sum_charge_ * RE * Gamma * Gamma * mol_density_ * density_correction;
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

    // v-dependence differs from bremsstrahlung
    return ((xi / 3)
               * (G / (Gamma * Gamma)
                   + 2 * ((x * x)  + (1 - x) * (1 - x)) * fi / Gamma))
        / (1 - (4. / 3) * x * (1 - x));
}
