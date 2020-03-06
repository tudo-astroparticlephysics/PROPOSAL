
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Ionization.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

/******************************************************************************
 *                               Ionization                                *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Ionization::Ionization(const ParticleDef& particle_def,
                       std::shared_ptr<const Medium> medium,
                       const EnergyCutSettings& cuts,
                       double multiplier)
    : Parametrization(particle_def, medium, cuts, multiplier)
{
}

Ionization::Ionization(const Ionization& ioniz)
    : Parametrization(ioniz)
{
}

Ionization::~Ionization() {}

// ------------------------------------------------------------------------- //
double Ionization::Delta(double beta, double gamma) {
    /* std::shared_ptr<const Medium> medium = this->GetMedium(); */
    double X;

    X = std::log(beta * gamma) / std::log(10);

    if (X < medium_->GetX0())
    {
        return medium_->GetD0() * std::pow(10, 2 * (X - medium_->GetX0()));
    } else if (X < medium_->GetX1())
    {
        return 2 * LOG10 * X + medium_->GetC() + medium_->GetA() * std::pow(medium_->GetX1() - X, medium_->GetM());
    } else
    {
        return 2 * LOG10 * X + medium_->GetC();
    }
}

// ------------------------------------------------------------------------- //
// Specific Parametrization
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits IonizBetheBlochRossi::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double mass_ration = ME / particle_def_.mass;

    double gamma = energy / particle_def_.mass;

    limits.vMin = (1.e-6 * medium_->GetI()) / energy;

    // PDG
    // eq. 33.4
    // v_{max} = \frac{1}{E} \frac{2 m_e \beta^2 \gamma^2}
    //          {1 + 2 \gamma \frac{m_e}{m_{particle} + (\frac{m_e}{m_{particle})^2 }
    limits.vMax = 2 * ME * (gamma * gamma - 1) / ((1 + 2 * gamma * mass_ration + mass_ration * mass_ration) * energy);
    limits.vMax = std::min(limits.vMax, 1. - particle_def_.mass / energy);

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

IonizBetheBlochRossi::IonizBetheBlochRossi(const ParticleDef& particle_def,
                                           std::shared_ptr<const Medium> medium,
                                           const EnergyCutSettings& cuts,
                                           double multiplier)
    : Ionization(particle_def, medium, cuts, multiplier)
{
}

IonizBetheBlochRossi::IonizBetheBlochRossi(const IonizBetheBlochRossi& ioniz)
    : Ionization(ioniz)
{
}

IonizBetheBlochRossi::~IonizBetheBlochRossi() {}

// ------------------------------------------------------------------------- //
// knonk-on electrons (delta rays)
// distribution of secondary electrons with kinetic energy = v*E
// PDG, Chin. Phys. C 40 (2016), 100001
// eq. 33.8
// ------------------------------------------------------------------------- //
double IonizBetheBlochRossi::DifferentialCrossSection(double energy, double v)
{
    double result;

    IntegralLimits limits = GetIntegralLimits(energy);

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum   = (energy - particle_def_.mass) * (energy + particle_def_.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta              = particle_momentum / energy;
    double gamma             = energy / particle_def_.mass;
    beta *= beta;

    // additional term for spin 1/2 particles
    // Rossi, 1952
    // High Enegy Particles
    // Prentice-Hall, Inc., Englewood Cliffs, N.J.
    // chapter 2, eq. 7
    double spin_1_2_contribution = v / (1 + 1 / gamma);
    spin_1_2_contribution *= 0.5 * spin_1_2_contribution;
    result = 1 - beta * (v / limits.vMax) + spin_1_2_contribution;
    result *= IONK * particle_def_.charge * particle_def_.charge * medium_->GetZA() / (2 * beta * energy * v * v);

    return medium_->GetMassDensity() * result * (1 + InelCorrection(energy, v));
    ;
}

// ------------------------------------------------------------------------- //
double IonizBetheBlochRossi::FunctionToDEdxIntegral(double energy, double variable)
{
    double result, aux;

    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);
    ParticleDef particle_def = this->GetParticleDef();
    /* std::shared_ptr<const Medium> medium     = this->GetMedium(); */

    // TODO(mario): Better way? Sat 2017/09/02

    // PDG eq. 33.10
    // with Spin 1/2 correction by Rossi
    double square_momentum   = (energy - particle_def.mass) * (energy + particle_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta              = particle_momentum / energy;
    double gamma             = energy / particle_def.mass;

    aux    = beta * gamma / (1.e-6 * medium_->GetI());
    result = std::log(limits.vUp * (2 * ME * energy)) + 2 * std::log(aux);
    aux    = limits.vUp / (2 * (1 + 1 / gamma));
    result += aux * aux;
    aux = beta * beta;
    result -= aux * (1 + limits.vUp / limits.vMax) + Delta(beta, gamma);

    if (result > 0)
    {
        result *= IONK * particle_def.charge * particle_def.charge * medium_->GetZA() / (2 * aux);
    } else
    {
        result = 0;
    }

    if(limits.vUp != limits.vMin){
        result *= medium_->GetMassDensity()/(limits.vUp - limits.vMin);
    }
    else{
        return 0;
    }

    return result / energy + variable * CrossSectionWithoutInelasticCorrection(energy, variable) * InelCorrection(energy, variable);
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
    IntegralLimits limits = GetIntegralLimits(energy);

    double gamma = energy / particle_def_.mass;

    a      = std::log(1 + 2 * v * energy / ME);
    b      = std::log((1 - v / limits.vMax) / (1 - v));
    c      = std::log((2 * gamma * (1 - v) * ME) / (particle_def_.mass * v));
    result = a * (2 * b + c) - b * b;

    return ALPHA / (2 * PI) * result;
}

// ------------------------------------------------------------------------- //
// CrossSection without inelastic correction
// needed for the dEdx Integral
// ------------------------------------------------------------------------- //

double IonizBetheBlochRossi::CrossSectionWithoutInelasticCorrection(double energy, double v)
{
    double result;

    IntegralLimits limits = GetIntegralLimits(energy);

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum   = (energy - particle_def_.mass) * (energy + particle_def_.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta              = particle_momentum / energy;
    double gamma             = energy / particle_def_.mass;
    beta *= beta;

    // additional term for spin 1/2 particles
    // Rossi, 1952
    // High Enegy Particles
    // Prentice-Hall, Inc., Englewood Cliffs, N.J.
    // chapter 2, eq. 7
    double spin_1_2_contribution = v / (1 + 1 / gamma);
    spin_1_2_contribution *= 0.5 * spin_1_2_contribution;
    result = 1 - beta * (v / limits.vMax) + spin_1_2_contribution;
    result *= IONK * particle_def_.charge * particle_def_.charge * medium_->GetZA() / (2 * beta * energy * v * v);

    return medium_->GetMassDensity() * result;
}


const std::string IonizBetheBlochRossi::name_ = "IonizBetheBlochRossi";

// ------------------------------------------------------------------------- //
// Ionization formula for positrons
// BetheBloch can't be used due to the ambiguity of the final state
// For high energy transfers, Moller or Bhabha Scattering is used (different for electron and positron)
// For low energies, we need to take a sum over the excitation probabilities of the atom (first order independent!)
// ------------------------------------------------------------------------- //

IonizBergerSeltzerBhabha::IonizBergerSeltzerBhabha(const ParticleDef& particle_def,
                                           std::shared_ptr<const Medium> medium,
                                           const EnergyCutSettings& cuts,
                                           double multiplier)
        : Ionization(particle_def, medium, cuts, multiplier)
{
}

IonizBergerSeltzerBhabha::IonizBergerSeltzerBhabha(const IonizBergerSeltzerBhabha& ioniz)
        : Ionization(ioniz)
{
}

IonizBergerSeltzerBhabha::~IonizBergerSeltzerBhabha() {}

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits IonizBergerSeltzerBhabha::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 0.;

    limits.vMax = 1. - ME / energy;

    if (limits.vMax < 0) {
        limits.vMax = 0;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}

double IonizBergerSeltzerBhabha::DifferentialCrossSection(double energy, double v)
{
    /*
     * Bhabha-Crosssection, taken from : "The EGS5 Code System",
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman, Scott and R. Nelson, Walter,
     * DOI: 10.2172/877459
     *
     * For high energy transfers, corresponding to high v and therefore stoachastic losses, the electrons of the
     * shell can be treated as free and we can use the Bhabha Scattering crosssection. Furthermore, we use the Bhabha
     * crosssection to estimate the dE2dx integral
     */

    double aux    = 0;

    double gamma = energy/ME;
    double epsilon = (v * energy)/(energy - ME);
    double betasquared = 1. - 1. / (gamma * gamma);
    double y = 1. / (gamma + 1.);
    double B1 = 2. - y * y;
    double B2 = (1. - 2. * y) * (3. + y * y);
    double B4 = std::pow(1. - 2. * y, 3.);
    double B3 = std::pow(1. - 2. * y, 2.) + B4;

    aux = 1. / (betasquared * epsilon * epsilon) - B1 / epsilon + B2 - B3 * epsilon + B4 * epsilon * epsilon;

    aux *= 1. / (gamma - 1.);
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI  * std::pow(RE, 2.) * medium_->GetMassDensity() * NA * medium_->GetZA() ;

    return std::max(aux, 0.);
}

// ------------------------------------------------------------------------- //
double IonizBergerSeltzerBhabha::FunctionToDEdxIntegral(double energy, double variable)
{
    /* Berger Seltzer Formula, taken from:
     * The EGS5 Code System,
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman, Scott and R. Nelson, Walter,
     * DOI: 10.2172/877459
     *
     * In general, especially for small energy transfers, one has to consider the excitation probabilities of the
     * atom. The Berger-Seltzer formula used here takes these information into account for small energy transfers
     * and uses the Bhabha-formula for high energy transfers, see:
     * "Positron-Electron Differences in Energy Loss and Multiple Scattering", F. Rohrlich and B.C. Carlson (1954)
     */

    (void)variable; // integral is calculated analytically here

    double result, aux;

    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);
    ParticleDef particle_def = this->GetParticleDef();
    std::shared_ptr<const Medium> medium     = this->GetMedium();

    double fplus; // (2.269)
    double bigDelta     = std::min(limits.vMax * energy / ME, limits.vUp * energy / ME); // (2.265)
    double gamma        = energy / ME; // (2.258)
    double betasquared  = 1. - 1. / (gamma * gamma); // (2.260)
    double tau          = gamma - 1.; // (2.262)
    double y            = 1. / (gamma + 1.); // (2.263)

    aux = tau + 2 * bigDelta - (3. * bigDelta * bigDelta * y / 2.) - (bigDelta - std::pow(bigDelta, 3.) / 3.) * y * y
            - (bigDelta * bigDelta / 2. - tau * std::pow(bigDelta, 3.) / 3. + std::pow(bigDelta, 4.) / 4.) * std::pow(y, 3.);

    aux *= betasquared / tau;
    fplus = std::log(tau * bigDelta) - aux;

    result = std::log( 2. * (tau + 2.) / (std::pow(1e-6 * medium_->GetI(), 2. ) / ME) );
    result += fplus;
    result -= Delta(std::sqrt(betasquared), gamma);

    result *= 2. * PI * RE * RE * ME / betasquared;

    result *= NA * medium_->GetZA() * medium_->GetMassDensity();

    return std::max(result, 0.);
}

const std::string IonizBergerSeltzerBhabha::name_ = "IonizBergerSeltzerBhabha";

// ------------------------------------------------------------------------- //
// Ionization formula for electrons
// BetheBloch can't be used due to the ambiguity of the final state
// For high energy transfers, Moller or Bhabha Scattering is used (different for electron and positron)
// For low energies, we need to take a sum over the excitation probabilities of the atom (first order independent!)
// ------------------------------------------------------------------------- //

IonizBergerSeltzerMoller::IonizBergerSeltzerMoller(const ParticleDef& particle_def,
                                                   std::shared_ptr<const Medium> medium,
                                                   const EnergyCutSettings& cuts,
                                                   double multiplier)
        : Ionization(particle_def, medium, cuts, multiplier)
{
}

IonizBergerSeltzerMoller::IonizBergerSeltzerMoller(const IonizBergerSeltzerMoller& ioniz)
        : Ionization(ioniz)
{
}

IonizBergerSeltzerMoller::~IonizBergerSeltzerMoller() {}

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits IonizBergerSeltzerMoller::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 0.;

    limits.vMax = 0.5 * (1. - ME / energy);

    if (limits.vMax < 0) {
        limits.vMax = 0;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}

double IonizBergerSeltzerMoller::DifferentialCrossSection(double energy, double v)
{

    /*
     * Moller-Crosssection, taken from : "The EGS5 Code System",
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman, Scott and R. Nelson, Walter,
     * DOI: 10.2172/877459
     * For high energy transfers, corresponding to high v and therefore stoachastic losses, the electrons of the
     * shell can be treated as free and we can use the Moller Scattering crosssection. Furthermore, we use the Moller
     * crosssection to estimate the dE2dx integral
     */

    double aux    = 0;

    double gamma = energy/ME;
    double epsilon = (v * energy)/(energy - ME);
    double betasquared = 1. - 1. / (gamma * gamma);

    aux = std::pow(gamma - 1., 2.) / (gamma * gamma) + 1. / epsilon * (1. / epsilon - (2. * gamma - 1.) / (gamma * gamma))
            + 1. / (1. - epsilon) * (1. / (1. - epsilon) - (2. * gamma - 1.) / (gamma * gamma));
    aux = aux / betasquared;

    aux *= 1. / (gamma - 1.);
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI  * std::pow(RE, 2.) * medium_->GetMassDensity() * NA * medium_->GetZA() ;

    return std::max(aux, 0.);
}

// ------------------------------------------------------------------------- //
double IonizBergerSeltzerMoller::FunctionToDEdxIntegral(double energy, double variable)
{
    /* Berger Seltzer Formula, taken from:
     * The EGS5 Code System,
     * Hirayama, Hideo and Namito, Yoshihito and Bielajew, Alex and Wilderman, Scott and R. Nelson, Walter,
     * DOI: 10.2172/877459
     *
     * In general, especially for small energy transfers, one has to consider the excitation probabilities of the
     * atom. The Berger-Seltzer formula used here takes these information into account for small energy transfers
     * and uses the Moller-formula for high energy transfers, see:
     * "Positron-Electron Differences in Energy Loss and Multiple Scattering", F. Rohrlich and B.C. Carlson (1954)
     */

    (void)variable; // integral is calculated analytically here

    double result, aux;

    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);
    ParticleDef particle_def = this->GetParticleDef();
    std::shared_ptr<const Medium> medium     = this->GetMedium();

    double fminus; // (2.268)
    double bigDelta     = std::min(limits.vMax * energy / ME, limits.vUp * energy / ME); // (2.265)
    double gamma        = energy / ME; // (2.258)
    double betasquared  = 1. - 1. / (gamma * gamma); // (2.260)
    double tau          = gamma - 1.; // (2.262)

    aux = bigDelta * bigDelta / 2. + (2. * tau + 1.) * std::log(1. - bigDelta / tau);
    aux = aux / (gamma * gamma);

    fminus = aux - 1. - betasquared + std::log((tau - bigDelta) * bigDelta) + tau / (tau - bigDelta);

    result = std::log( 2. * (tau + 2.) / (std::pow(1e-6 * medium->GetI(), 2. ) / ME) );
    result += fminus;
    result -= Delta(std::sqrt(betasquared), gamma);

    result *= 2. * PI * RE * RE * ME / betasquared;

    result *= NA * medium->GetZA() * medium_->GetMassDensity();

    return std::max(result, 0.);
}

const std::string IonizBergerSeltzerMoller::name_ = "IonizBergerSeltzerMoller";
