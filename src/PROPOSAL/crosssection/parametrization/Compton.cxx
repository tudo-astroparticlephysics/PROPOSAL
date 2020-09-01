
#include <cmath>
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;
using std::make_tuple;

crosssection::Compton::Compton()
    : Parametrization(InteractionType::Compton, "compton")
{
}

double crosssection::Compton::GetLowerEnergyLim(const ParticleDef&) const noexcept
{
    return ME;
}

tuple<double, double> crosssection::Compton::GetKinematicLimits(
    const ParticleDef&, const Component&, double energy) const noexcept
{
    assert(energy > 0);
    auto vmax = 1. - 1. / (1. + 2. * energy / ME);
    return make_tuple(0., vmax);
}

double crosssection::ComptonKleinNishina::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy, double v) const
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

namespace PROPOSAL {
template <>
double integrate_dndx(Integral& integral, crosssection::Compton& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max, double rnd)
{
    auto t_min = std::log(1. - v_max);
    auto t_max = std::log(1. - v_min);
    auto dNdx = [&param, &p_def, &comp, energy](double t) {
        return exp(t)
            * param.FunctionToDNdxIntegral(p_def, comp, energy, 1 - exp(t));
    };
    return integral.IntegrateWithRandomRatio(t_min, t_max, dNdx, 4, rnd);
}

template <>
double integrate_dedx(Integral& integral, crosssection::Compton& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max)
{
    auto t_min = std::log(1. - v_max);
    auto t_max = std::log(1. - v_min);
    auto dEdx = [&param, &p_def, &comp, energy](double t) {
        return exp(t)
            * param.FunctionToDEdxIntegral(p_def, comp, energy, 1 - exp(t));
    };
    return integral.Integrate(t_min, t_max, dEdx, 2);
}

template <>
double integrate_de2dx(Integral& integral, crosssection::Compton& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max)
{
    auto t_min = std::log(1. - v_max);
    auto t_max = std::log(1. - v_min);
    auto dE2dx = [&param, &p_def, &comp, energy](double t) {
        return exp(t)
            * param.FunctionToDE2dxIntegral(p_def, comp, energy, 1 - exp(t));
    };
    return integral.Integrate(t_min, t_max, dE2dx, 2);
}

/* template <> */
/* double calculate_upper_lim_dndx(Integral& integral, crosssection::Compton& param, */
/*     const ParticleDef& p_def, const Component& comp, double energy, */
/*     double v_min, double v_max, double rnd) */
/* { */
/*     auto t_min = std::log(1. - v_min); */
/*     auto t_max = std::log(1. - v_max); */
/*     auto dNdx = [&param, &p_def, &comp, energy](double t) { */
/*         return exp(t) */
/*             * param.FunctionToDNdxIntegral(p_def, comp, energy, 1 - exp(t)); */
/*     }; */
/*     integral.IntegrateWithRandomRatio(t_min, t_max, dNdx, 4, rnd); */
/*     return integral.GetUpperLimit(); */
/* } */
} // namespace PROPOSAL
