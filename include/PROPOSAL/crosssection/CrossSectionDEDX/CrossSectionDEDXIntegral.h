#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {

namespace crosssection {
    class IonizBergerSeltzerBhabha;
    class IonizBergerSeltzerMoller;
}

namespace detail {
    using dedx_integral_t = std::function<double(Integral&, double)>;

    template <typename T1>
    dedx_integral_t define_dedx_integral(T1 param, ParticleDef const& p_def,
        Component const& comp, EnergyCutSettings const& cut)
    {
        return [param, p_def, comp, cut](Integral& i, double E) {
            auto lim = param.GetKinematicLimits(p_def, comp, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [&param, &p_def, &comp, E](double v) {
                return param.FunctionToDEdxIntegral(p_def, comp, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 2);
        };
    }

    template <typename T1>
    dedx_integral_t define_dedx_integral(T1 param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const& cut)
    {
        return [param, p_def, medium, cut](Integral& i, double E) {
            auto lim = param.GetKinematicLimits(p_def, medium, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [&param, &p_def, &medium, E](double v) {
                return param.FunctionToDEdxIntegral(p_def, medium, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 4);
        };
    }

    template <>
    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&);

    template <>
    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&);
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
