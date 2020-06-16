#pragma once

#include "PROPOSAL/EnergyCutSettings.h"

#include <memory>

using std::shared_ptr;

namespace PROPOSAL {
using v_trafo_t = std::function<double(double, double, double)>;

class CrossSectionDNDX {
protected:
    std::function<double(double, double)> function_to_dndx;
    std::function<tuple<double, double>(double)> kinematic_limits;
    std::function<tuple<double, double>(double)> integral_limits;

    size_t hash_cross_section;

public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDX(const Param& _param, Particle _p_def, Target _target,
        shared_ptr<EnergyCutSettings> _cut)
        : function_to_dndx([_param, _p_def, _target](double energy, double v) {
            return _param.FunctionToDNdxIntegral(_p_def, _target, energy, v);
        })
        , kinematic_limits([_param, _p_def, _target](double energy) {
            return _param.GetKinematicLimits(_p_def, _target, energy);
        })
        , integral_limits([this, _cut](double energy) {
            return GetIntegralLimits(_cut.get(), energy);
        })
    {
        hash_combine(hash_cross_section, _param.GetHash(), _p_def.GetHash(),
            _target.GetHash());
    }

    virtual double Calculate(double, double, v_trafo_t = nullptr) = 0;
    virtual double GetUpperLim(double, double, v_trafo_t = nullptr) = 0;

    enum { MIN, MAX };
    tuple<double, double> GetIntegralLimits(EnergyCutSettings*, double);
    tuple<double, double> GetIntegralLimits(std::nullptr_t, double);
};

} // namespace PROPOSAL
