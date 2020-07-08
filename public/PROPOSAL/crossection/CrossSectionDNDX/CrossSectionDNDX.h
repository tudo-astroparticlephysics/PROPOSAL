#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"

#include <memory>

using std::function;
using std::make_tuple;
using std::shared_ptr;

namespace PROPOSAL {
using v_trafo_t = std::function<double(double, double, double)>;

class CrossSectionDNDX {
protected:
    function<double(Integral&, double, double, double, double)> dndx_integral;
    function<tuple<double, double>(double)> kinematic_limits;
    shared_ptr<const EnergyCutSettings> cut;
    size_t hash_cross_section;

public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDX(Param _param, Particle _particle, Target _target,
        shared_ptr<const EnergyCutSettings> _cut)
        : dndx_integral(
              [_param, _particle, _target](Integral& integral, double energy,
                  double v_min, double v_max, double rate) {
                  return integrate_dndx(integral, _param, _particle, _target,
                      energy, v_min, v_max, rate);
              })
        , kinematic_limits([_param, _particle, _target](double energy) {
            return _param.GetKinematicLimits(_particle, _target, energy);
        })
        , cut(_cut)
    {
        hash_cross_section = 0;
        hash_combine(hash_cross_section, _param.GetHash(), _particle.GetHash(),
            _target.GetHash());
    }

    virtual ~CrossSectionDNDX() = default;

    inline double Calculate(double energy)
    {
        return Calculate(energy, get<MAX>(GetIntegrationLimits(energy)));
    }
    virtual double Calculate(double energy, double v, v_trafo_t = nullptr) = 0;
    virtual double GetUpperLimit(
        double energy, double rate, v_trafo_t = nullptr)
        = 0;

    enum { MIN, MAX };
    inline tuple<double, double> GetIntegrationLimits(double energy)
    {
        auto kin_lim = kinematic_limits(energy);
        if (cut) {
            return make_tuple(cut->GetCut(kin_lim, energy),
                              get<crosssection::Parametrization::V_MAX>(kin_lim));
        } else {
            return make_tuple(get<crosssection::Parametrization::V_MIN>(kin_lim),
                              get<crosssection::Parametrization::V_MAX>(kin_lim));
        }

    }
};
} // namespace PROPOSAL
