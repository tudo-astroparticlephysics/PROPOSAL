
#include <cmath>
#include <stdexcept>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

using std::get;
using std::logic_error;
using std::make_tuple;
using std::tuple;
using namespace PROPOSAL;

crosssection::Ionization::Ionization(const EnergyCutSettings& cuts)
    : cuts_(cuts)
{
    hash_combine(hash, cuts_.GetHash());
}

double crosssection::Ionization::GetLowerEnergyLim(
    const ParticleDef& p_def) const noexcept
{
    return p_def.mass;
}

double crosssection::Ionization::Delta(
    const Medium& medium, double beta, double gamma) const
{
    auto X = std::log(beta * gamma) / std::log(10);

    if (X < medium.GetX0()) {
        return medium.GetD0() * std::pow(10, 2 * (X - medium.GetX0()));
    } else if (X < medium.GetX1()) {
        return 2 * LOG10 * X + medium.GetC()
            + medium.GetA() * std::pow(medium.GetX1() - X, medium.GetM());
    } else {
        return 2 * LOG10 * X + medium.GetC();
    }
}

crosssection::KinematicLimits
crosssection::IonizBetheBlochRossi::GetKinematicLimits(
    const ParticleDef& p_def, const Medium& medium, double energy) const
{
    auto mass_ration = ME / p_def.mass;
    auto gamma = energy / p_def.mass;
    auto kin_lim = KinematicLimits();
    kin_lim.v_min = (1.e-6 * medium.GetI()) / energy;
    // PDG eq. 33.4
    // v_{max} = \frac{1}{E} \frac{2 m_e \beta^2 \gamma^2}
    //          {1 + 2 \gamma \frac{m_e}{m_{particle} +
    //          (\frac{m_e}{m_{particle})^2 }
    kin_lim.v_max = 2 * ME * (gamma * gamma - 1)
        / ((1 + 2 * gamma * mass_ration + mass_ration * mass_ration) * energy);
    kin_lim.v_max = std::min(kin_lim.v_max, 1. - p_def.mass / energy);
    if (kin_lim.v_max < kin_lim.v_min)
        kin_lim.v_max = kin_lim.v_min;
    return kin_lim;
}

crosssection::IonizBetheBlochRossi::IonizBetheBlochRossi(
    const EnergyCutSettings& cuts)
    : crosssection::Ionization(cuts)
{
    hash_combine(hash, std::string("bethe_bloch_rossi"));
}

std::unique_ptr<crosssection::Parametrization<Medium>>
crosssection::IonizBetheBlochRossi::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// knonk-on electrons (delta rays)
// distribution of secondary electrons with kinetic energy = v*E
// PDG, Chin. Phys. C 40 (2016), 100001
// eq. 33.8
// ------------------------------------------------------------------------- //
double crosssection::IonizBetheBlochRossi::DifferentialCrossSection(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double v) const
{

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = (energy - p_def.mass) * (energy + p_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta = particle_momentum / energy;
    double gamma = energy / p_def.mass;
    beta *= beta;

    // additional term for spin 1/2 particles
    // Rossi, 1952
    // High Enegy Particles
    // Prentice-Hall, Inc., Englewood Cliffs, N.J.
    // chapter 2, eq. 7
    double spin_1_2_contribution = v / (1 + 1 / gamma);
    spin_1_2_contribution *= 0.5 * spin_1_2_contribution;
    auto result = 1
        - beta * (v / GetKinematicLimits(p_def, medium, energy).v_max)
        + spin_1_2_contribution;
    result *= IONK * p_def.charge * p_def.charge
        * calculate_proton_massnumber_fraction(medium.GetComponents())
        / (2 * beta * energy * v * v);

    return result * (1 + InelCorrection(p_def, medium, energy, v));
}

// ------------------------------------------------------------------------- //
double crosssection::IonizBetheBlochRossi::FunctionToDEdxIntegral(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double variable) const
{
    return variable
      * CrossSectionWithoutInelasticCorrection(p_def, medium, energy, variable)
      * InelCorrection(p_def, medium, energy, variable);
}

double crosssection::IonizBetheBlochRossi::IonizationLoss(
        const ParticleDef& p_def, const Medium& medium, double energy) const {
    double result, aux;

    auto limits = GetKinematicLimits(p_def, medium, energy);

    // TODO(mario): Better way? Sat 2017/09/02
    // PDG eq. 33.10
    // with Spin 1/2 correction by Rossi
    double square_momentum = (energy - p_def.mass) * (energy + p_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta = particle_momentum / energy;
    double gamma = energy / p_def.mass;
    double v_up;

    v_up = limits.v_max;
    v_up = std::min(v_up, cuts_.GetCut(energy));

    aux = beta * gamma / (1.e-6 * medium.GetI());
    result = std::log(v_up * (2 * ME * energy)) + 2 * std::log(aux);
    aux = v_up / (2 * (1 + 1 / gamma));
    result += aux * aux;
    aux = beta * beta;
    result -= aux * (1 + v_up / limits.v_max) + Delta(medium, beta, gamma);

    if (result > 0) {
        result *= IONK * p_def.charge * p_def.charge
                  * calculate_proton_massnumber_fraction(medium.GetComponents())
                  / (2 * aux);
    } else {
        result = 0;
    }

    if (v_up == limits.v_min)
        return 0;
    return result;
}

// ------------------------------------------------------------------------- //
// Bremststrahlung when scattering at atomic electrons
// and the atomic electrons emit the Bremsstrahlung photon
// because of the v^{-2} dependency, it is treated together with Ionization
// Kelner Kokoulin Petrukhin
// Phys. Atom. Nucl. 60 (1997), 657
// eq. 30
// \Delta \frac{d \sigma}{d v} = \frac{d \sigma}{d v}_{I_0}
//     \frac{\alpha}{2 \pi} \cdot
//        (\log(1 + \frac{2vE}{m_e})
//        (2 \log(\frac{1 - \frac{v}{v_{max}}}{1 - v}))
//        \log(\frac{2 \gamma (1 - v) m_e}{m_{particle}v})
// ------------------------------------------------------------------------- //
double crosssection::IonizBetheBlochRossi::InelCorrection(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double v) const
{
    double gamma = energy / p_def.mass;
    auto a = std::log(1 + 2 * v * energy / ME);
    auto b = std::log(
        (1 - v / GetKinematicLimits(p_def, medium, energy).v_max) / (1 - v));
    auto c = std::log((2 * gamma * (1 - v) * ME) / (p_def.mass * v));
    auto result = a * (2 * b + c) - b * b;

    return ALPHA / (2 * PI) * result;
}

// ------------------------------------------------------------------------- //
// CrossSection without inelastic correction
// needed for the dEdx Integral
// ------------------------------------------------------------------------- //

double
crosssection::IonizBetheBlochRossi::CrossSectionWithoutInelasticCorrection(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double v) const
{
    double result;

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = (energy - p_def.mass) * (energy + p_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta = particle_momentum / energy;
    double gamma = energy / p_def.mass;
    beta *= beta;

    // additional term for spin 1/2 particles
    // Rossi, 1952
    // High Enegy Particles
    // Prentice-Hall, Inc., Englewood Cliffs, N.J.
    // chapter 2, eq. 7
    double spin_1_2_contribution = v / (1 + 1 / gamma);
    spin_1_2_contribution *= 0.5 * spin_1_2_contribution;
    result = 1 - beta * (v / GetKinematicLimits(p_def, medium, energy).v_max)
        + spin_1_2_contribution;
    result *= IONK * p_def.charge * p_def.charge
        * calculate_proton_massnumber_fraction(medium.GetComponents())
        / (2 * beta * energy * v * v);

    return result;
}

// ------------------------------------------------------------------------- //
// Ionization formula for positrons
// BetheBloch can't be used due to the ambiguity of the final state
// For high energy transfers, Moller or Bhabha Scattering is used (different for
// electron and positron) For low energies, we need to take a sum over the
// excitation probabilities of the atom (first order independent!)
// ------------------------------------------------------------------------- //

crosssection::IonizBergerSeltzerBhabha::IonizBergerSeltzerBhabha(
    const EnergyCutSettings& cuts)
    : crosssection::Ionization(cuts)
{
    hash_combine(hash, std::string("berger_seltzer_bhabha"));
}

std::unique_ptr<crosssection::Parametrization<Medium>>
crosssection::IonizBergerSeltzerBhabha::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

crosssection::KinematicLimits
crosssection::IonizBergerSeltzerBhabha::GetKinematicLimits(
    const ParticleDef& p_def, const Medium&, double energy) const
{
    auto kin_lim = KinematicLimits();
    kin_lim.v_min = 0.;
    kin_lim.v_max = 1. - p_def.mass / energy;

    if (kin_lim.v_max < 0)
        kin_lim.v_max = 0;

    return kin_lim;
}

double crosssection::IonizBergerSeltzerBhabha::DifferentialCrossSection(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double v) const
{
    /*
     * Bhabha-Crosssection, taken from : "The EGS5 Code System",
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman,
     * Scott and R. Nelson, Walter, DOI: 10.2172/877459
     *
     * For high energy transfers, corresponding to high v and therefore
     * stoachastic losses, the electrons of the shell can be treated as free and
     * we can use the Bhabha Scattering crosssection. Furthermore, we use the
     * Bhabha crosssection to estimate the dE2dx integral
     */

    double aux = 0;

    double gamma = energy / p_def.mass;
    double epsilon = (v * energy) / (energy - p_def.mass);
    double betasquared = 1. - 1. / (gamma * gamma);
    double y = 1. / (gamma + 1.);
    double B1 = 2. - y * y;
    double B2 = (1. - 2. * y) * (3. + y * y);
    double B4 = std::pow(1. - 2. * y, 3.);
    double B3 = std::pow(1. - 2. * y, 2.) + B4;

    aux = 1. / (betasquared * epsilon * epsilon) - B1 / epsilon + B2
        - B3 * epsilon + B4 * epsilon * epsilon;

    aux *= 1. / (gamma - 1.);
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI * std::pow(RE, 2.) * NA
        * calculate_proton_massnumber_fraction(medium.GetComponents());

    return std::max(aux, 0.);
}

// ------------------------------------------------------------------------- //
double crosssection::IonizBergerSeltzerBhabha::FunctionToDEdxIntegral(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double variable) const
{
    /* Berger Seltzer Formula, taken from:
     * The EGS5 Code System,
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman,
     * Scott and R. Nelson, Walter, DOI: 10.2172/877459
     *
     * In general, especially for small energy transfers, one has to consider
     * the excitation probabilities of the atom. The Berger-Seltzer formula used
     * here takes these information into account for small energy transfers and
     * uses the Bhabha-formula for high energy transfers, see:
     * "Positron-Electron Differences in Energy Loss and Multiple Scattering",
     * F. Rohrlich and B.C. Carlson (1954)
     */

    (void)variable; // integral is calculated analytically here

    double result, aux;
    double fplus; // (2.269)
    double v_up;
    auto limits = GetKinematicLimits(p_def, medium, energy);

    v_up = limits.v_max;
    v_up = std::min(cuts_.GetCut(energy), v_up);

    double bigDelta = std::min(limits.v_max * energy / ME,
        v_up * energy / ME);                        // (2.265)
    double gamma = energy / ME;                     // (2.258)
    double betasquared = 1. - 1. / (gamma * gamma); // (2.260)
    double tau = gamma - 1.;                        // (2.262)
    double y = 1. / (gamma + 1.);                   // (2.263)

    aux = tau + 2 * bigDelta - (3. * bigDelta * bigDelta * y / 2.)
        - (bigDelta - std::pow(bigDelta, 3.) / 3.) * y * y
        - (bigDelta * bigDelta / 2. - tau * std::pow(bigDelta, 3.) / 3.
              + std::pow(bigDelta, 4.) / 4.)
            * std::pow(y, 3.);

    aux *= betasquared / tau;
    fplus = std::log(tau * bigDelta) - aux;

    result
        = std::log(2. * (tau + 2.) / (std::pow(1e-6 * medium.GetI(), 2.) / ME));
    result += fplus;
    result -= Delta(medium, std::sqrt(betasquared), gamma);

    result *= 2. * PI * RE * RE * ME / betasquared;

    result *= NA * calculate_proton_massnumber_fraction(medium.GetComponents());

    return std::max(result, 0.);
}

// ------------------------------------------------------------------------- //
// Ionization formula for electrons
// BetheBloch can't be used due to the ambiguity of the final state
// For high energy transfers, Moller or Bhabha Scattering is used (different for
// electron and positron) For low energies, we need to take a sum over the
// excitation probabilities of the atom (first order independent!)
// ------------------------------------------------------------------------- //

crosssection::IonizBergerSeltzerMoller::IonizBergerSeltzerMoller(
    const EnergyCutSettings& cuts)
    : crosssection::Ionization(cuts)
{
    hash_combine(hash, std::string("berger_seltzer_moller"));
}

std::unique_ptr<crosssection::Parametrization<Medium>>
crosssection::IonizBergerSeltzerMoller::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

crosssection::KinematicLimits
crosssection::IonizBergerSeltzerMoller::GetKinematicLimits(
    const ParticleDef& p_def, const Medium&, double energy) const
{
    auto kin_lim = KinematicLimits();
    kin_lim.v_min = 0.;
    kin_lim.v_max = 0.5 * (1. - p_def.mass / energy);

    if (kin_lim.v_max < 0)
        kin_lim.v_max = 0;

    return kin_lim;
}

double crosssection::IonizBergerSeltzerMoller::DifferentialCrossSection(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double v) const
{
    /*
     * Moller-Crosssection, taken from : "The EGS5 Code System",
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman,
     * Scott and R. Nelson, Walter, DOI: 10.2172/877459 For high energy
     * transfers, corresponding to high v and therefore stoachastic losses, the
     * electrons of the shell can be treated as free and we can use the Moller
     * Scattering crosssection. Furthermore, we use the Moller crosssection to
     * estimate the dE2dx integral
     */

    double aux = 0;

    double gamma = energy / p_def.mass;
    double epsilon = (v * energy) / (energy - p_def.mass);
    double betasquared = 1. - 1. / (gamma * gamma);

    aux = std::pow(gamma - 1., 2.) / (gamma * gamma)
        + 1. / epsilon * (1. / epsilon - (2. * gamma - 1.) / (gamma * gamma))
        + 1. / (1. - epsilon)
            * (1. / (1. - epsilon) - (2. * gamma - 1.) / (gamma * gamma));
    aux = aux / betasquared;

    aux *= 1. / (gamma - 1.);
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI * std::pow(RE, 2.) * NA
        * calculate_proton_massnumber_fraction(medium.GetComponents());

    return std::max(aux, 0.);
}

// ------------------------------------------------------------------------- //
double crosssection::IonizBergerSeltzerMoller::FunctionToDEdxIntegral(
    const ParticleDef& p_def, const Medium& medium, double energy,
    double variable) const
{
    /* Berger Seltzer Formula, taken from:
     * The EGS5 Code System,
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman,
     * Scott and R. Nelson, Walter, DOI: 10.2172/877459
     *
     * In general, especially for small energy transfers, one has to consider
     * the excitation probabilities of the atom. The Berger-Seltzer formula used
     * here takes these information into account for small energy transfers and
     * uses the Moller-formula for high energy transfers, see:
     * "Positron-Electron Differences in Energy Loss and Multiple Scattering",
     * F. Rohrlich and B.C. Carlson (1954)
     */

    (void)variable; // integral is calculated analytically here

    double result, aux;

    auto limits = GetKinematicLimits(p_def, medium, energy);

    double fminus; // (2.268)
    double v_up;

    v_up = limits.v_max;
    v_up = std::min(cuts_.GetCut(energy), v_up);

    double bigDelta = std::min(limits.v_max * energy / ME,
        v_up * energy / ME);                        // (2.265)
    double gamma = energy / ME;                     // (2.258)
    double betasquared = 1. - 1. / (gamma * gamma); // (2.260)
    double tau = gamma - 1.;                        // (2.262)

    aux = bigDelta * bigDelta / 2.
        + (2. * tau + 1.) * std::log(1. - bigDelta / tau);
    aux = aux / (gamma * gamma);

    fminus = aux - 1. - betasquared + std::log((tau - bigDelta) * bigDelta)
        + tau / (tau - bigDelta);

    result
        = std::log(2. * (tau + 2.) / (std::pow(1e-6 * medium.GetI(), 2.) / ME));
    result += fminus;
    result -= Delta(medium, std::sqrt(betasquared), gamma);

    result *= 2. * PI * RE * RE * ME / betasquared;

    result *= NA * calculate_proton_massnumber_fraction(medium.GetComponents());

    return std::max(result, 0.);
}
