#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL {
class CrossSectionDEDXIntegral : public CrossSectionDEDX {
    Integral integral;
    std::function<double(Integral&, double)> dedx_integral;

    /* template <typename... Args> inline auto define_integral(Args... args) */
    /* { */
    /*     return [args...](Integral& i, double E, double v_min, double v_max) {
     */
    /*         return integrate_dedx(i, args..., E, v_min, v_max); */
    /*     }; */
    /* } */

    template <typename T1, typename T2>
    inline auto define_integral(T1 param, ParticleDef const& p_def, T2 comp,
        EnergyCutSettings const& cut)
    {
        return [param, p_def, comp, cut](Integral& i, double E) {
            auto lim = param.GetKinematicLimits(p_def, comp, E);
            auto v_cut = cut.GetCut(lim, E);
            return integrate_dedx(i, param, p_def, comp, E,
                std::get<crosssection::Parametrization::V_MIN>(lim), v_cut);
        };
    }

public:
    template <typename... Args>
    CrossSectionDEDXIntegral(Args... args)
        : dedx_integral(define_integral(std::forward<Args>(args)...))
    {
    }

    double Calculate(double E) override { return dedx_integral(integral, E); }
};
} // namespace PROPOSAL
