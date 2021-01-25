
#include <cmath>
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

double crosssection::Compton::GetLowerEnergyLim(
    const ParticleDef&) const noexcept
{
    return ME;
}

std::tuple<double, double> crosssection::Compton::GetKinematicLimits(
    const ParticleDef&, const Component&, double energy) const noexcept
{
    assert(energy > 0);
    auto vmax = 1. - 1. / (1. + 2. * energy / ME);
    return std::make_tuple(0., vmax);
}

double crosssection::ComptonKleinNishina::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    // Adapted from "THE EGS5 CODE SYSTEM" by Hideo Harayama and Yoshihito
    // Namito SLAC Report number SLAC-R-730, KEK Report number 2005-8 Equation
    // (2.170) adapted for PROPOSAL Based on O. Klein and Y. Nishina "Über die
    // Streuung von Strahlung durch freie Electronen nach der Neuen
    // Relativistischen Quantum Dynamic von Dirac", Zeitschrift für Physik, 52,
    // 853-868.
    assert(energy >= p_def.mass);
    (void)p_def;
    auto kp = energy / ME;
    auto C1 = std::pow(kp, -2.);
    auto C2 = 1. - 2. * (1. + kp) * C1;
    auto C3 = (1 + 2. * kp) * C1;
    auto epsilon = 1. - v;
    auto aux = (C1 / epsilon + C2) / epsilon + C3 + epsilon;
    aux /= energy; // we loose a factor E due to variable transformation from k'
                   // to v
    aux *= PI * std::pow(RE, 2.) * ME;
    return NA / comp.GetAtomicNum() * comp.GetNucCharge() * aux;
}
