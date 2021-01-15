#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL {

namespace detail {
    template <typename T1, typename T2, typename T3>
    auto define_dndx_integral(T1 param, T2 p_def, T3 target)
    {
        return [param, p_def, target](Integral& i, double E, double v_min,
                   double v_max, double rate) {
            auto dNdx = [&param, &p_def, &target, E](double v) {
                return param.DifferentialCrossSection(p_def, target, E, v);
            };
            return i.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rate);
        };
    }

    template <typename T1, typename T2, typename T3>
    auto define_dndx_upper_lim(T1 param, T2 p_def, T3 target)
    {
        return [param, p_def, target](Integral& i, double E, double v_min,
                   double v_max, double rnd) {
            auto dNdx = [&param, &p_def, &target, E](double v) {
                return param.DifferentialCrossSection(p_def, target, E, v);
            };
            i.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rnd);
            return i.GetUpperLimit();
        };
    }
}

class CrossSectionDNDXIntegral : public CrossSectionDNDX {
    Integral integral;
    std::function<double(Integral&, double, double, double, double)>
        dndx_integral;
    std::function<double(Integral&, double, double, double, double)>
        dndx_upper_limit;

public:
    template <typename T1, typename T2, typename T3, typename... Args>
    CrossSectionDNDXIntegral(T1 _param, T2 _particle, T3 _target, Args... args)
        : CrossSectionDNDX(_param, _particle, _target, args...)
        , dndx_integral(
              detail::define_dndx_integral(_param, _particle, _target))
        , dndx_upper_limit(
              detail::define_dndx_upper_lim(_param, _particle, _target))
    {
    }

    double Calculate(double energy) override;

    double Calculate(double energy, double v) override;

    double GetUpperLimit(double energy, double rate) override;
};
} // namespace PROPOSAL
