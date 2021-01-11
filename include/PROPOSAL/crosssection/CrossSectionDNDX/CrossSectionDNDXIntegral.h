#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL {
class CrossSectionDNDXIntegral : public CrossSectionDNDX {
    Integral integral;
    std::function<double(Integral&, double, double, double, double)>
        dndx_integral;
    std::function<double(Integral&, double, double, double, double)>
        dndx_upper_limit;

    template <typename... Args> inline auto define_integral(Args... args)
    {
        return [args...](Integral& i, double E, double v_min, double v_max,
                   double rate) {
            return integrate_dndx(i, args..., E, v_min, v_max, rate);
        };
    }

    template <typename... Args> inline auto define_upper_lim(Args... args)
    {
        return [args...](Integral& i, double E, double v_min, double v_max,
                   double rnd) {
            return calculate_upper_lim_dndx(i, args..., E, v_min, v_max, rnd);
        };
    }

public:
    template <typename T1, typename T2, typename T3, typename... Args>
    CrossSectionDNDXIntegral(T1 _param, T2 _particle, T3 _target, Args... args)
        : CrossSectionDNDX(_param, _particle, _target, args...)
        , dndx_integral(define_integral(_param, _particle, _target))
        , dndx_upper_limit(define_upper_lim(_param, _particle, _target))
    {
    }

    double Calculate(double energy) override;

    double Calculate(double energy, double v) override;

    double GetUpperLimit(double energy, double rate) override;
};
} // namespace PROPOSAL
