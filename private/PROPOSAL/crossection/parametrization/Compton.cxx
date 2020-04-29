
#include <cmath>
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

Compton::Compton(const ParticleDef& p_def, const component_list& comp)
    : Parametrization("compton", p_def, comp, ME)
{
}

Parametrization::KinematicLimits Compton::GetKinematicLimits(double energy)
{
    assert(energy > 0);
    auto vmax = 1. - 1. / (1. + 2. * energy / ME);

    return KinematicLimits(0., vmax);
}

ComptonKleinNishina::ComptonKleinNishina(
    const ParticleDef& p_def, const component_list& comp)
    : Compton(p_def, comp)
{
}

double ComptonKleinNishina::DifferentialCrossSection(double energy, double v)
{
    // Adapted from "THE EGS5 CODE SYSTEM" by Hideo Harayama and Yoshihito
    // Namito SLAC Report number SLAC-R-730, KEK Report number 2005-8 Equation
    // (2.170) adapted for PROPOSAL Based on O. Klein and Y. Nishina "Über die
    // Streuung von Strahlung durch freie Electronen nach der Neuen
    // Relativistischen Quantum Dynamic von Dirac", Zeitschrift für Physik, 52,
    // 853-868.

    assert(energy >= particle_mass_);

    auto kp = energy / ME;
    auto C1 = std::pow(kp, -2.);
    auto C2 = 1. - 2. * (1. + kp) * C1;
    auto C3 = (1 + 2. * kp) * C1;
    auto epsilon = 1. - v;

    auto aux = (C1 / epsilon + C2) / epsilon + C3 + epsilon;

    aux /= energy; // we loose a factor E due to variable transformation from k'
                   // to v

    aux *= PI * std::pow(RE, 2.) * ME;
    return current_component_.GetAtomInMolecule()
        * current_component_.GetNucCharge() * aux;
}
