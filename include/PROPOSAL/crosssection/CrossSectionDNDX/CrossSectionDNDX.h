#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"

#include <memory>

namespace PROPOSAL {
/* using v_trafo_t = std::function<double(double, double, double)>; */

class CrossSectionDNDX {
    template <typename T1, typename T2, typename T3>
    inline auto define_kinematic_limits(T1 param, T2 particle, T3 target)
    {
        return [param, particle, target](double E) {
            return param.GetKinematicLimits(particle, target, E);
        };
    }

    std::function<std::tuple<double, double>(double)> kinematic_limits;
    std::shared_ptr<const EnergyCutSettings> cut;
    double lower_energy_lim;

public:
    template <typename T1, typename T2, typename T3>
    CrossSectionDNDX(T1 _param, T2 _particle, T3 _target,
        std::shared_ptr<const EnergyCutSettings> _cut)
        : kinematic_limits(define_kinematic_limits(_param, _particle, _target))
        , cut(_cut)
        , lower_energy_lim(_param.GetLowerEnergyLim(_particle))
    {
    }

    virtual ~CrossSectionDNDX() = default;

    virtual double Calculate(double energy) = 0;
    virtual double Calculate(double energy, double v) = 0;
    virtual double GetUpperLimit(double energy, double rate) = 0;

    enum { MIN, MAX };
    inline std::array<double, 2> GetIntegrationLimits(double energy)
    {
        auto kin_lim = kinematic_limits(energy);
        if (cut)
            return std::array<double, 2> { cut->GetCut(kin_lim, energy),
                std::get<crosssection::Parametrization::V_MAX>(kin_lim) };
        return std::array<double, 2> {
            std::get<crosssection::Parametrization::V_MIN>(kin_lim),
            std::get<crosssection::Parametrization::V_MAX>(kin_lim)
        };
    }

    auto GetLowerEnergyLim() { return lower_energy_lim; }
};
} // namespace PROPOSAL
