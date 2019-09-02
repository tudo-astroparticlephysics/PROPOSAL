
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
                       const Medium& medium,
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
Parametrization::IntegralLimits Ionization::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double mass_ration = ME / particle_def_.mass;
    ;
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

// ------------------------------------------------------------------------- //
// Specific Parametrization
// ------------------------------------------------------------------------- //

IonizBetheBlochRossi::IonizBetheBlochRossi(const ParticleDef& particle_def,
                                           const Medium& medium,
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
    const Medium& medium     = this->GetMedium();

    // TODO(mario): Better way? Sat 2017/09/02

    // PDG eq. 33.10
    // with Spin 1/2 correction by Rossi
    double square_momentum   = (energy - particle_def.mass) * (energy + particle_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta              = particle_momentum / energy;
    double gamma             = energy / particle_def.mass;

    aux    = beta * gamma / (1.e-6 * medium.GetI());
    result = std::log(limits.vUp * (2 * ME * energy)) + 2 * std::log(aux);
    aux    = limits.vUp / (2 * (1 + 1 / gamma));
    result += aux * aux;
    aux = beta * beta;
    result -= aux * (1 + limits.vUp / limits.vMax) + Delta(beta, gamma);

    if (result > 0)
    {
        result *= IONK * particle_def.charge * particle_def.charge * medium.GetZA() / (2 * aux);
    } else
    {
        result = 0;
    }

    if(limits.vUp != limits.vMin){
        result *= medium.GetMassDensity()/(limits.vUp - limits.vMin);
    }
    else{
        return 0;
    }

    return result / energy + variable * CrossSectionWithoutInelasticCorrection(energy, variable) * InelCorrection(energy, variable);
}

// ------------------------------------------------------------------------- //
double IonizBetheBlochRossi::Delta(double beta, double gamma) {
    const Medium& medium = this->GetMedium();
    double X;

    X = std::log(beta * gamma) / std::log(10);

    if (X < medium.GetX0())
    {
        return medium.GetD0() * std::pow(10, 2 * (X - medium.GetX0()));
    } else if (X < medium.GetX1())
    {
        return 2 * LOG10 * X + medium.GetC() + medium.GetA() * std::pow(medium.GetX1() - X, medium.GetM());
    } else
    {
        return 2 * LOG10 * X + medium.GetC();
    }
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
