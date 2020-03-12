
#include <functional>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

/******************************************************************************
 *                               Compton                                      *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Compton::Compton(const ParticleDef& particle_def,
                               std::shared_ptr<const Medium> medium,
                               const EnergyCutSettings& cuts,
                               double multiplier)
        : Parametrization(particle_def, medium, cuts, multiplier)
{
}

Compton::Compton(const Compton& compton)
        : Parametrization(compton)
{
}

Compton::~Compton() {}

bool Compton::compare(const Parametrization& parametrization) const
{
    //const Compton* compton = static_cast<const Compton*>(&parametrization);
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Compton::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 0.;

    // Limits from kinematics

    double aux = 1. - 1. / (1. + 2. * energy / ME);

    if(aux > 1){
        limits.vMax = 1;
    }
    else{
        limits.vMax = aux;
    }

    if (limits.vMax < 0)
    {
        limits.vMax = 0;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}


// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

size_t Compton::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed);

    return seed;
}

/******************************************************************************
 *                          Specifc Parametrizations                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

ComptonKleinNishina::ComptonKleinNishina(const ParticleDef& particle_def,
                                               std::shared_ptr<const Medium> medium,
                                               const EnergyCutSettings& cuts,
                                               double multiplier)
        : Compton(particle_def, medium, cuts, multiplier)
{
}

ComptonKleinNishina::ComptonKleinNishina(const ComptonKleinNishina& compton)
        : Compton(compton)
{
}

ComptonKleinNishina::~ComptonKleinNishina() {
}

bool ComptonKleinNishina::compare(const Parametrization& parametrization) const
{
    //const ComptonKleinNishina* compton = static_cast<const ComptonKleinNishina*>(&parametrization);
    return Compton::compare(parametrization);
}

const std::string ComptonKleinNishina::name_ = "ComptonKleinNishina";

double ComptonKleinNishina::DifferentialCrossSection(double energy, double v)
{
    // Adapted from "THE EGS5 CODE SYSTEM" by Hideo Harayama and Yoshihito Namito
    // SLAC Report number SLAC-R-730, KEK Report number 2005-8
    // Equation (2.170) adapted for PROPOSAL
    // Based on O. Klein and Y. Nishina "Über die Streuung von Strahlung durch freie
    // Electronen nach der Neuen Relativistischen Quantum Dynamic von Dirac", Zeitschrift für Physik, 52, 853-868.

    double aux    = 0;

    double C1, C2, C3, epsilon;
    double kp = energy / ME;

    C1 = std::pow(kp, -2.);
    C2 = 1. - 2. * (1. + kp) * C1;
    C3 = (1 + 2. * kp) * C1;
    epsilon = 1. - v;

    aux = ( C1 / epsilon + C2 ) / epsilon + C3 + epsilon;

    aux = aux / energy; // we loose a factor E due to variable transformation from k' to v

    aux *= PI * std::pow(RE, 2.) * ME;
    return medium_->GetMolDensity() * components_[component_index_].GetAtomInMolecule() * components_[component_index_].GetNucCharge() * aux;
}
