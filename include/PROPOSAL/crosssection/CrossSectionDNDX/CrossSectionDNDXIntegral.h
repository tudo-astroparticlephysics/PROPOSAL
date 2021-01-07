#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include <memory>
#include <unordered_map>

using PROPOSAL::Components::Component;

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
    template <typename T1, typename T2, typename T3>
    CrossSectionDNDXIntegral(T1 _param, T2 _particle, T3 _target,
        std::shared_ptr<const EnergyCutSettings> _cut)
        : CrossSectionDNDX(_param, _particle, _target, _cut)
        , dndx_integral(define_integral(_param, _particle, _target))
        , dndx_upper_limit(define_upper_lim(_param, _particle, _target))
    {
    }

    double Calculate(double energy) override
    {
        return Calculate(energy, std::get<MAX>(GetIntegrationLimits(energy)));
    }

    double Calculate(double energy, double v) override
    {
        auto integral_lim = GetIntegrationLimits(energy);
        if (std::get<MIN>(integral_lim) < v)
            return dndx_integral(
                integral, energy, std::get<MIN>(integral_lim), v, 0);
        return 0;
    }

    double GetUpperLimit(double energy, double rate) override
    {
        auto integral_lim = GetIntegrationLimits(energy);
        auto v = dndx_upper_limit(integral, energy, std::get<MIN>(integral_lim),
            std::get<MAX>(integral_lim), -rate);
        return v;
    }
};
} // namespace PROPOSAL
