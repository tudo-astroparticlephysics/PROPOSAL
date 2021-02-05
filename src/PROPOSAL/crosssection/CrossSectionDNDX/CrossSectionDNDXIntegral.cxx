#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

namespace PROPOSAL {
namespace detail {

    template <typename T> using param_t = crosssection::Parametrization<T>;

    template <typename Target>
    dndx_integrand_t _define_dndx_integral(
        param_t<Target> const& param, ParticleDef p, Target t)
    {
        return [ptr = std::shared_ptr<param_t<Target>>(param.clone()), p, t](
                   double E, double v_min, double v_max, double rate) {
            Integral i;
            auto dNdx = [param_ptr = ptr.get(), &p, &t, E](double v) {
                return param_ptr->DifferentialCrossSection(p, t, E, v);
            };
            return i.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rate);
        };
    }

    dndx_integrand_t define_dndx_integral(
        param_t<Medium> const& param, ParticleDef p, Medium m)
    {
        return _define_dndx_integral(param, p, m);
    }

    dndx_integrand_t define_dndx_integral(
        param_t<Component> const& param, ParticleDef p, Component c)
    {
        return _define_dndx_integral(param, p, c);
    }

    template <typename Target>
    dndx_upper_lim_t _define_dndx_upper_lim(
        param_t<Target> const& param, ParticleDef p, Target t)
    {
        return [ptr = std::shared_ptr<param_t<Target>>(param.clone()), p, t](
                   double E, double v_min, double v_max, double rnd) {
            Integral i;
            auto dNdx = [param_ptr = ptr.get(), &p, &t, E](double v) {
                return param_ptr->DifferentialCrossSection(p, t, E, v);
            };
            i.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rnd);
            return i.GetUpperLimit();
        };
    }

    dndx_upper_lim_t define_dndx_upper_lim(
        crosssection::Parametrization<Medium> const& param, ParticleDef p,
        Medium m)
    {
        return _define_dndx_upper_lim(param, p, m);
    }

    dndx_upper_lim_t define_dndx_upper_lim(
        crosssection::Parametrization<Component> const& param, ParticleDef p,
        Component c)
    {
        return _define_dndx_upper_lim(param, p, c);
    }

} // namespace PROPOSAL
} // namespace detail

double CrossSectionDNDXIntegral::Calculate(double energy)
{
    return Calculate(energy, std::get<MAX>(GetIntegrationLimits(energy)));
}

double CrossSectionDNDXIntegral::Calculate(double energy, double v)
{
    auto lim = GetIntegrationLimits(energy);
    if (std::get<MIN>(lim) < v)
        return dndx_integral(energy, std::get<MIN>(lim), v, 0);
    return 0;
}

double CrossSectionDNDXIntegral::GetUpperLimit(double energy, double rate)
{
    auto lim = GetIntegrationLimits(energy);
    auto v = dndx_upper_limit(
        energy, std::get<MIN>(lim), std::get<MAX>(lim), -rate);
    return v;
}
