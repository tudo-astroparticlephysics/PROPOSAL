#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
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
        param_t<Target> const& param, ParticleDef const& p, Target const& t)
    {
        return [ptr = std::shared_ptr<param_t<Target>>(param.clone()), p, t](
                   double E, double v_min, double v_max) {
            Integral i;
            auto dNdx = [param_ptr = ptr.get(), &p, &t, E](double v) {
                return param_ptr->DifferentialCrossSection(p, t, E, v);
            };
            return i.Integrate(v_min, v_max, dNdx, 4);
        };
    }

    dndx_integrand_t define_dndx_integral(
        param_t<Medium> const& param, ParticleDef const& p, Medium const& m)
    {
        return _define_dndx_integral(param, p, m);
    }

    dndx_integrand_t define_dndx_integral(param_t<Component> const& param,
        ParticleDef const& p, Component const& c)
    {
        return _define_dndx_integral(param, p, c);
    }

    dndx_integrand_t define_dndx_integral(
        crosssection::ComptonKleinNishina const& param, ParticleDef const& p,
        Component const& c)
    {
        using param_t = crosssection::Parametrization<Component>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, c](double E, double v_min, double v_max) {
            Integral i;
            auto dNdx = [ptr = param_ptr.get(), &p, &c, E](double v) {
                double t = std::log(1. - v);
                return ptr->DifferentialCrossSection(p, c, E, t);
            };
            return i.Integrate(v_min, v_max, dNdx, 2);
        };
    }

    template <typename Target>
    dndx_upper_lim_t _define_dndx_upper_lim(
        param_t<Target> const& param, ParticleDef const& p, Target const& t)
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
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p, Medium const& m)
    {
        return _define_dndx_upper_lim(param, p, m);
    }

    dndx_upper_lim_t define_dndx_upper_lim(
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p, Component const& c)
    {
        return _define_dndx_upper_lim(param, p, c);
    }

    dndx_upper_lim_t define_dndx_upper_lim(
        crosssection::ComptonKleinNishina const& param, ParticleDef const& p,
        Component const& c)
    {
        using param_t = crosssection::Parametrization<Component>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, c](
                   double E, double v_min, double v_max, double rate) {
            Integral i;
            auto dNdx = [ptr = param_ptr.get(), &p, &c, E](double v) {
                double t = std::log(1. - v);
                return ptr->DifferentialCrossSection(p, c, E, t);
            };
            i.IntegrateWithRandomRatio(v_min, v_max, dNdx, 3, rate);
            return i.GetUpperLimit();
        };
    }
} // namespace PROPOSAL
} // namespace detail

double CrossSectionDNDXIntegral::Calculate(double energy)
{
    auto lim = GetIntegrationLimits(energy);
    return Calculate(energy, lim.max);
}

double CrossSectionDNDXIntegral::Calculate(double energy, double v)
{
    auto lim = GetIntegrationLimits(energy);
    if (lim.min < v)
        return dndx_integral(energy, lim.min, v);
    return 0;
}

double CrossSectionDNDXIntegral::GetUpperLimit(double energy, double rate)
{
    auto lim = GetIntegrationLimits(energy);
    auto v = dndx_upper_limit(energy, lim.min, lim.max, -rate);
    return v;
}
