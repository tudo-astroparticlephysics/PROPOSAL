#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <memory>

using namespace PROPOSAL;

namespace PROPOSAL {
namespace detail {

    template <typename Param, typename Target>
    dedx_integral_t _define_dedx_integral(Param const& param,
        ParticleDef const& p, Target const& t, EnergyCutSettings const& cut)
    {
        auto param_ptr = std::shared_ptr<Param>(param.clone());
        return [param_ptr, p, t, cut](double E) {
            auto i = Integral();
            auto lim = param_ptr->GetKinematicLimits(p, t, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [_param_ptr = param_ptr.get(), &p, &t, E](double v) {
                return _param_ptr->FunctionToDEdxIntegral(p, t, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 2);
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p, Component const& c, EnergyCutSettings const& cut)
    {
        return _define_dedx_integral(param, p, c, cut);
    }

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut)
    {
        return _define_dedx_integral(param, p, m, cut);
    }

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const&)
    {
        using param_t = crosssection::Parametrization<Medium>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, m](double E) {
            return param_ptr->FunctionToDEdxIntegral(p, m, E, 0.) / E;
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const&)
    {
        using param_t = crosssection::Parametrization<Medium>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, m](double E) {
            return param_ptr->FunctionToDEdxIntegral(p, m, E, 0.) / E;
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::ComptonKleinNishina const& param, ParticleDef const& p,
        Component const& c, EnergyCutSettings const& cut)
    {
        using param_t = crosssection::Parametrization<Component>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, c, cut](double E) {
            // Integrate with the substitution t = ln(1-v) to avoid
            // numerical problems
            auto i = Integral();
            auto lim = param_ptr->GetKinematicLimits(p, c, E);
            auto v_cut = cut.GetCut(lim, E);
            double t_min = std::log(1. - v_cut);
            double t_max = std::log(1. - lim.v_min);
            auto dEdx = [ptr = param_ptr.get(), p, c, E](double t) {
                return std::exp(t)
                    * ptr->FunctionToDEdxIntegral(p, c, E, 1 - std::exp(t));
            };
            return i.Integrate(t_min, t_max, dEdx, 2);
        };
    }
}
}

double CrossSectionDEDXIntegral::Calculate(double E) const
{
    return dedx_integral(E) * E;
}
