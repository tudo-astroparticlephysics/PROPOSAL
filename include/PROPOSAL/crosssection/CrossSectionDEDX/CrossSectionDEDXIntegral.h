#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {

namespace detail {
    template <typename T1>
    auto define_dedx_integral(T1 param, ParticleDef const& p_def,
        Component const& comp, EnergyCutSettings const& cut)
    {
        return [param, p_def, comp, cut](Integral& i, double E) {
            auto lim = param.GetKinematicLimits(p_def, comp, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [&param, &p_def, &comp, E](double v) {
                return param.FunctionToDEdxIntegral(p_def, comp, E, v);
            };
            return i.Integrate(
                std::get<crosssection::Parametrization::V_MIN>(lim), v_cut,
                dEdx, 2);
        };
    }

    template <typename T1>
    auto define_dedx_integral(T1 param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const& cut)
    {
        return [param, p_def, medium, cut](Integral& i, double E) {
            auto physical_lim = param.GetKinematicLimits(p_def, medium, E);
            auto v_cut = cut.GetCut(physical_lim, E);
            auto dEdx = [&param, &p_def, &medium, E](double v) {
                return param.FunctionToDEdxIntegral(p_def, medium, E, v);
            };
            return i.Integrate(
                std::get<crosssection::Parametrization::V_MIN>(physical_lim),
                v_cut, dEdx, 4);
        };
    }

}

class CrossSectionDEDXIntegral : public CrossSectionDEDX {
    Integral integral;
    std::function<double(Integral&, double)> dedx_integral;

public:
    template <typename... Args>
    CrossSectionDEDXIntegral(Args... args)
        : CrossSectionDEDX(args...)
        , dedx_integral(detail::define_dedx_integral(args...))
    {
    }

    double Calculate(double E) final { return dedx_integral(integral, E) * E; }
};

} // namespace PROPOSAL
