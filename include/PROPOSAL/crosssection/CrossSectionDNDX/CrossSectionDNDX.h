#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"

#include <memory>

namespace PROPOSAL {
using v_trafo_t = std::function<double(double, double, double)>;

class CrossSectionDNDX {
protected:
    std::function<double(Integral&, double, double, double, double)> dndx_integral;
    std::function<double(Integral&, double, double, double, double)> dndx_upper_limit;
    std::function<std::tuple<double, double>(double)> kinematic_limits;
    std::shared_ptr<const EnergyCutSettings> cut;
    size_t hash_cross_section;

public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDX(Param _param, Particle _particle, Target _target,
        std::shared_ptr<const EnergyCutSettings> _cut)
        : dndx_integral(
              [_param, _particle, _target](Integral& integral, double energy,
                  double v_min, double v_max, double rate) mutable {
                  using param_t = typename std::decay<Param>::type;
                  using base_param_ref_t = typename std::add_lvalue_reference<
                          typename param_t::base_param_t>::type;
                  return integrate_dndx(
                          integral, reinterpret_cast<base_param_ref_t>(_param),
                          _particle, _target, energy, v_min, v_max, rate);
              }),
          dndx_upper_limit(
                  [_param, _particle, _target](Integral& integral, double energy,
                                               double v_min, double v_max, double rnd) mutable {
                      using param_t = typename std::decay<Param>::type;
                      using base_param_ref_t = typename std::add_lvalue_reference<
                              typename param_t::base_param_t>::type;
                      return calculate_upper_lim_dndx(
                              integral, reinterpret_cast<base_param_ref_t>(_param),
                              _particle, _target, energy, v_min, v_max, rnd);
                  })
        , kinematic_limits([_param, _particle, _target](double energy) {
            return _param.GetKinematicLimits(_particle, _target, energy);
        })
        , cut(_cut)
    {
        hash_cross_section = 0;
        hash_combine(hash_cross_section, _param.GetHash(), _particle.mass,
                     std::abs(_particle.charge), _target.GetHash());
        if (_param.interaction_type == InteractionType::WeakInt)
            hash_combine(hash_cross_section, _particle.charge);
        if (cut)
            hash_combine(hash_cross_section, cut->GetHash());
    }

    virtual ~CrossSectionDNDX() = default;

    virtual double Calculate(double energy) = 0;
    virtual double Calculate(double energy, double v) = 0;
    virtual double GetUpperLimit(double energy, double rate) = 0;

    enum { MIN, MAX };
    inline std::tuple<double, double> GetIntegrationLimits(double energy)
    {
        auto kin_lim = kinematic_limits(energy);
        if (cut) {
            return std::make_tuple(cut->GetCut(kin_lim, energy),
                              std::get<crosssection::Parametrization::V_MAX>(kin_lim));
        } else {
            return std::make_tuple(std::get<crosssection::Parametrization::V_MIN>(kin_lim),
                              std::get<crosssection::Parametrization::V_MAX>(kin_lim));
        }

    }
};
} // namespace PROPOSAL
