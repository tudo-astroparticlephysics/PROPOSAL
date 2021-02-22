#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXIntegral.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

namespace PROPOSAL {
namespace detail {
    template <typename Param, typename Target>
    std::function<double(double)> _define_de2dx_integral(
        Param const& param, ParticleDef p, Target t, EnergyCutSettings cut)
    {
        auto param_ptr = std::shared_ptr<Param>(param.clone());
        return [param_ptr, p, t, cut](double E) {
            auto i = Integral();
            auto lim = param_ptr->GetKinematicLimits(p, t, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dE2dx = [_param_ptr = param_ptr.get(), &p, &t, E](double v) {
                return _param_ptr->FunctionToDE2dxIntegral(p, t, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dE2dx, 2);
        };
    }

    std::function<double(double)> define_de2dx_integral(
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut)
    {
        return _define_de2dx_integral(param, p, m, cut);
    }

    std::function<double(double)> define_de2dx_integral(
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p, Component const& c, EnergyCutSettings const& cut)
    {
        return _define_de2dx_integral(param, p, c, cut);
    }
} // namespace detail
} // namespace PROPOSAL

double CrossSectionDE2DXIntegral::Calculate(double E) const
{
    return de2dx_integral(E) * E * E;
}
