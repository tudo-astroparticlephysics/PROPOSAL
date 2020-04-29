
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

Ionization::Ionization(const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cuts)
    : Parametrization("ionization", p_def, medium.GetComponents(), p_def.mass)
    , cuts_(cuts)
    , X0_(medium.GetX0())
    , X1_(medium.GetX1())
    , D0_(medium.GetD0())
    , m_(medium.GetM())
    , C_(medium.GetC())
    , a_(medium.GetA())
    , I_(medium.GetI())
    , ZA_(calculate_proton_massnumber_fraction(medium.GetComponents()))
{
}

double Ionization::Delta(double beta, double gamma)
{
    auto X = std::log(beta * gamma) / std::log(10);

    if (X < X0_) {
        return D0_ * std::pow(10, 2 * (X - X0_));
    } else if (X < X1_) {
        return 2 * LOG10 * X + C_ + a_ * std::pow(X1_ - X, m_);
    } else {
        return 2 * LOG10 * X + C_;
    }
}

Parametrization::KinematicLimits IonizBetheBlochRossi::GetKinematicLimits(
    double energy)
{
    auto mass_ration = ME / particle_mass_;
    auto gamma = energy / particle_mass_;

    auto v_min = (1.e-6 * I_) / energy;

    // PDG eq. 33.4
    // v_{max} = \frac{1}{E} \frac{2 m_e \beta^2 \gamma^2}
    //          {1 + 2 \gamma \frac{m_e}{m_{particle} +
    //          (\frac{m_e}{m_{particle})^2 }
    auto v_max = 2 * ME * (gamma * gamma - 1)
        / ((1 + 2 * gamma * mass_ration + mass_ration * mass_ration) * energy);
    v_max = std::min(v_max, 1. - particle_mass_ / energy);

    if (v_max < v_min)
        v_max = v_min;

    return KinematicLimits(v_min, v_max);
}

IonizBetheBlochRossi::IonizBetheBlochRossi(const ParticleDef& particle_def,
    const Medium& medium, const EnergyCutSettings& cuts)
    : Ionization(particle_def, medium, cuts)
{
}

// ------------------------------------------------------------------------- //
// knonk-on electrons (delta rays)
// distribution of secondary electrons with kinetic energy = v*E
// PDG, Chin. Phys. C 40 (2016), 100001
// eq. 33.8
// ------------------------------------------------------------------------- //
double IonizBetheBlochRossi::DifferentialCrossSection(double energy, double v)
{
    double result;

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum
        = (energy - particle_mass_) * (energy + particle_mass_);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta = particle_momentum / energy;
    double gamma = energy / particle_mass_;
    beta *= beta;

    // additional term for spin 1/2 particles
    // Rossi, 1952
    // High Enegy Particles
    // Prentice-Hall, Inc., Englewood Cliffs, N.J.
    // chapter 2, eq. 7
    double spin_1_2_contribution = v / (1 + 1 / gamma);
    spin_1_2_contribution *= 0.5 * spin_1_2_contribution;
    result = 1 - beta * (v / GetKinematicLimits(energy).vMax)
        + spin_1_2_contribution;
    result *= IONK * particle_charge_ * particle_charge_ * ZA_
        / (2 * beta * energy * v * v);

    return result * (1 + InelCorrection(energy, v));
    ;
}

// ------------------------------------------------------------------------- //
double IonizBetheBlochRossi::FunctionToDEdxIntegral(
    double energy, double variable)
{
    double result, aux;

    auto limits = GetKinematicLimits(energy);

    // TODO(mario): Better way? Sat 2017/09/02
    // PDG eq. 33.10
    // with Spin 1/2 correction by Rossi
    double square_momentum
        = (energy - particle_mass_) * (energy + particle_mass_);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta = particle_momentum / energy;
    double gamma = energy / particle_mass_;
    double v_up;

    v_up = limits.vMax;
    v_up = std::min(v_up, cuts_.GetCut(energy));

    aux = beta * gamma / (1.e-6 * I_);
    result = std::log(v_up * (2 * ME * energy)) + 2 * std::log(aux);
    aux = v_up / (2 * (1 + 1 / gamma));
    result += aux * aux;
    aux = beta * beta;
    result -= aux * (1 + v_up / limits.vMax) + Delta(beta, gamma);

    if (result > 0) {
        result *= IONK * particle_charge_ * particle_charge_ * ZA_ / (2 * aux);
    } else {
        result = 0;
    }

    if (v_up == limits.vMin)
        return 0;

    result /= (v_up - limits.vMin);

    return result / energy
        + variable * CrossSectionWithoutInelasticCorrection(energy, variable)
        * InelCorrection(energy, variable);
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
double IonizBetheBlochRossi::InelCorrection(double energy, double v)
{
    double result, a, b, c;

    double gamma = energy / particle_mass_;

    a = std::log(1 + 2 * v * energy / ME);
    b = std::log((1 - v / GetKinematicLimits(energy).vMax) / (1 - v));
    c = std::log((2 * gamma * (1 - v) * ME) / (particle_mass_ * v));
    result = a * (2 * b + c) - b * b;

    return ALPHA / (2 * PI) * result;
}

// ------------------------------------------------------------------------- //
// CrossSection without inelastic correction
// needed for the dEdx Integral
// ------------------------------------------------------------------------- //

double IonizBetheBlochRossi::CrossSectionWithoutInelasticCorrection(
    double energy, double v)
{
    double result;

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum
        = (energy - particle_mass_) * (energy + particle_mass_);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta = particle_momentum / energy;
    double gamma = energy / particle_mass_;
    beta *= beta;

    // additional term for spin 1/2 particles
    // Rossi, 1952
    // High Enegy Particles
    // Prentice-Hall, Inc., Englewood Cliffs, N.J.
    // chapter 2, eq. 7
    double spin_1_2_contribution = v / (1 + 1 / gamma);
    spin_1_2_contribution *= 0.5 * spin_1_2_contribution;
    result = 1 - beta * (v / GetKinematicLimits(energy).vMax)
        + spin_1_2_contribution;
    result *= IONK * particle_charge_ * particle_charge_ * ZA_
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

IonizBergerSeltzerBhabha::IonizBergerSeltzerBhabha(
    const ParticleDef& particle_def, const Medium& medium,
    const EnergyCutSettings& cuts)
    : Ionization(particle_def, medium, cuts)
{
}

Parametrization::KinematicLimits IonizBergerSeltzerBhabha::GetKinematicLimits(
    double energy)
{
    auto v_min = 0.;
    auto v_max = 1. - ME / energy;

    if (v_max < 0)
        v_max = 0;

    return KinematicLimits(v_min, v_max);
}

double IonizBergerSeltzerBhabha::DifferentialCrossSection(
    double energy, double v)
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

    double gamma = energy / ME;
    double epsilon = (v * energy) / (energy - ME);
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
    aux *= 2. * PI * std::pow(RE, 2.) * NA * ZA_;

    return std::max(aux, 0.);
}

// ------------------------------------------------------------------------- //
double IonizBergerSeltzerBhabha::FunctionToDEdxIntegral(
    double energy, double variable)
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
    auto limits = GetKinematicLimits(energy);

    v_up = limits.vMax;
    v_up = std::min(cuts_.GetCut(energy), v_up);

    double bigDelta
        = std::min(limits.vMax * energy / ME, v_up * energy / ME); // (2.265)
    double gamma = energy / ME;                                    // (2.258)
    double betasquared = 1. - 1. / (gamma * gamma);                // (2.260)
    double tau = gamma - 1.;                                       // (2.262)
    double y = 1. / (gamma + 1.);                                  // (2.263)

    aux = tau + 2 * bigDelta - (3. * bigDelta * bigDelta * y / 2.)
        - (bigDelta - std::pow(bigDelta, 3.) / 3.) * y * y
        - (bigDelta * bigDelta / 2. - tau * std::pow(bigDelta, 3.) / 3.
              + std::pow(bigDelta, 4.) / 4.)
            * std::pow(y, 3.);

    aux *= betasquared / tau;
    fplus = std::log(tau * bigDelta) - aux;

    result = std::log(2. * (tau + 2.) / (std::pow(1e-6 * I_, 2.) / ME));
    result += fplus;
    result -= Delta(std::sqrt(betasquared), gamma);

    result *= 2. * PI * RE * RE * ME / betasquared;

    result *= NA * ZA_;

    return std::max(result, 0.);
}

// ------------------------------------------------------------------------- //
// Ionization formula for electrons
// BetheBloch can't be used due to the ambiguity of the final state
// For high energy transfers, Moller or Bhabha Scattering is used (different for
// electron and positron) For low energies, we need to take a sum over the
// excitation probabilities of the atom (first order independent!)
// ------------------------------------------------------------------------- //

IonizBergerSeltzerMoller::IonizBergerSeltzerMoller(
    const ParticleDef& particle_def, const Medium& medium,
    const EnergyCutSettings& cuts)
    : Ionization(particle_def, medium, cuts)
{
}

Parametrization::KinematicLimits IonizBergerSeltzerMoller::GetKinematicLimits(
    double energy)
{
    auto v_min = 0.;
    auto v_max = 0.5 * (1. - ME / energy);

    if (v_max < 0)
        v_max = 0;

    return KinematicLimits(v_min, v_max);
}

double IonizBergerSeltzerMoller::DifferentialCrossSection(
    double energy, double v)
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

    double gamma = energy / ME;
    double epsilon = (v * energy) / (energy - ME);
    double betasquared = 1. - 1. / (gamma * gamma);

    aux = std::pow(gamma - 1., 2.) / (gamma * gamma)
        + 1. / epsilon * (1. / epsilon - (2. * gamma - 1.) / (gamma * gamma))
        + 1. / (1. - epsilon)
            * (1. / (1. - epsilon) - (2. * gamma - 1.) / (gamma * gamma));
    aux = aux / betasquared;

    aux *= 1. / (gamma - 1.);
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI * std::pow(RE, 2.) * NA * ZA_;

    return std::max(aux, 0.);
}

// ------------------------------------------------------------------------- //
double IonizBergerSeltzerMoller::FunctionToDEdxIntegral(
    double energy, double variable)
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

    auto limits = GetKinematicLimits(energy);

    double fminus; // (2.268)
    double v_up;

    v_up = limits.vMax;
    v_up = std::min(cuts_.GetCut(energy), v_up);

    double bigDelta
        = std::min(limits.vMax * energy / ME, v_up * energy / ME); // (2.265)
    double gamma = energy / ME;                                    // (2.258)
    double betasquared = 1. - 1. / (gamma * gamma);                // (2.260)
    double tau = gamma - 1.;                                       // (2.262)

    aux = bigDelta * bigDelta / 2.
        + (2. * tau + 1.) * std::log(1. - bigDelta / tau);
    aux = aux / (gamma * gamma);

    fminus = aux - 1. - betasquared + std::log((tau - bigDelta) * bigDelta)
        + tau / (tau - bigDelta);

    result = std::log(2. * (tau + 2.) / (std::pow(1e-6 * I_, 2.) / ME));
    result += fminus;
    result -= Delta(std::sqrt(betasquared), gamma);

    result *= 2. * PI * RE * RE * ME / betasquared;

    result *= NA * ZA_;

    return std::max(result, 0.);
}
