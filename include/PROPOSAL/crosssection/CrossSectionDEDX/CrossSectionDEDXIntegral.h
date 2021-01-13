#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {
class CrossSectionDEDXIntegral : public CrossSectionDEDX {
    Integral integral;
    std::function<double(Integral&, double)> dedx_integral;

    template <typename T1>
    auto define_integral(T1 param, ParticleDef const& p_def,
        Component const& comp, EnergyCutSettings const& cut)
    {
        return [param, p_def, comp, cut](Integral& i, double E) {
            auto lim = param.GetKinematicLimits(p_def, comp, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [&param, &p_def, &comp, E](double v) {
                return param.FunctionToDEdxIntegral(p_def, comp, E, v);
            };
            return i.Integrate(
                       std::get<crosssection::Parametrization::V_MIN>(lim),
                       v_cut, dEdx, 2)
                * E;
        };
    }

    template <typename T1>
    auto define_integral(T1 param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const& cut)
    {
        return [param, p_def, medium, cut](Integral& i, double E) {
            auto physical_lim = param.GetKinematicLimits(p_def, medium, E);
            auto v_cut = cut.GetCut(physical_lim, E);
            auto dEdx = [&param, &p_def, &medium, E](double v) {
                return E * param.FunctionToDEdxIntegral(p_def, medium, E, v);
            };
            return i.Integrate(
                std::get<crosssection::Parametrization::V_MIN>(physical_lim),
                v_cut, dEdx, 4);
        };
    }

    template <typename T1, typename T2>
    inline auto define_integral(T1, ParticleDef const&, T2, std::nullptr_t)
    {
        return [](Integral&, double) { return 0.; };
    }

public:
    template <typename... Args>
    CrossSectionDEDXIntegral(Args... args)
        : dedx_integral(define_integral(std::forward<Args>(args)...))
    {
    }

    double Calculate(double E) override { return dedx_integral(integral, E); }
};

/* namespace crosssection { */
/*  class IonizBergerSeltzerMoller; */
/*  class IonizBergerSeltzerBhabha; */
/* } */

template <>
inline auto CrossSectionDEDXIntegral::define_integral(
    crosssection::IonizBergerSeltzerBhabha param, ParticleDef const& p_def,
    Medium const& medium, EnergyCutSettings const& cut)

{
    return [param, p_def, medium](Integral&, double E) {
        return param.FunctionToDEdxIntegral(p_def, medium, E, 0.);
    };
}

template <>
inline auto CrossSectionDEDXIntegral::define_integral(
    crosssection::IonizBergerSeltzerMoller param, ParticleDef const& p_def,
    Medium const& medium, EnergyCutSettings const& cut)
{
    return [param, p_def, medium](Integral&, double E) {
        return param.FunctionToDEdxIntegral(p_def, medium, E, 0.);
    };
}

} // namespace PROPOSAL
