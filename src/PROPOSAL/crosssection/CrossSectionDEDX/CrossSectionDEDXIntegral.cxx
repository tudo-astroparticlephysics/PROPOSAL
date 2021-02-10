#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {
namespace detail {

    template <typename Param, typename Target>
    dedx_integral_t _define_dedx_integral(Param const& param,
        ParticleDef const& p, Target const& t, EnergyCutSettings const& cut)
    {
        auto param_ptr = std::shared_ptr<Param>(param.clone());
        return [param_ptr, p, t, cut](Integral& i, double E) {
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

    dedx_integral_t define_dedx_integral(crosssection::IonizBetheBlochRossi param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut)
    {
        using param_t = crosssection::Parametrization<Medium>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, m, cut](Integral& i, double E) {
            auto lim = param_ptr->GetKinematicLimits(p, m, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [_param_ptr = param_ptr.get(), &p, &m, E](double v) {
                return _param_ptr->FunctionToDEdxIntegral(p, m, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 4);
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&)
    {
        using param_t = crosssection::Parametrization<Medium>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p_def, medium](Integral&, double E) {
            return param_ptr->FunctionToDEdxIntegral(p_def, medium, E, 0.) / E;
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&)
    {
        using param_t = crosssection::Parametrization<Medium>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p_def, medium](Integral&, double E) {
            return param_ptr->FunctionToDEdxIntegral(p_def, medium, E, 0.) / E;
        };
    }
}
double CrossSectionDEDXIntegral::Calculate(double E) const
{
    auto integral = Integral();
    return dedx_integral(integral, E) * E;
}
}
