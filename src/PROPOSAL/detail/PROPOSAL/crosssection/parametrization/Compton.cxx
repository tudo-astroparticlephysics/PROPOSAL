
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

double crosssection::Compton::GetLowerEnergyLim(
    const ParticleDef&) const noexcept
{
    return 2. * ME;
}

crosssection::KinematicLimits crosssection::Compton::GetKinematicLimits(
    ParticleDef const&, Component const&, double energy) const
{
    assert(energy > 0);
    auto kin_lim = KinematicLimits();
    kin_lim.v_min = 0;
    kin_lim.v_max = 1. - 1. / (1. + 2. * energy / ME);
    return kin_lim;
}

std::unique_ptr<crosssection::Parametrization<Component>>
crosssection::ComptonKleinNishina::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
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
