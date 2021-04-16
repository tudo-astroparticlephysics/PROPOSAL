#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relative_loss(double v_cut, double v_max, double v);
double retransform_relative_loss(double v_cut, double v_max, double v);
double transform_loss_log(double v_cut, double v_max, double v);
double retransform_loss_log(double v_cut, double v_max, double v);

template <typename Param> struct transform_loss {
    static std::function<double(double, double, double)> func;
};

template <typename Param> struct retransform_loss {
    static std::function<double(double, double, double)> func;
};

template <typename Param>
std::function<double(double, double, double)> transform_loss<Param>::func
    = [](double v_cut, double v_max, double v) {
          return transform_relative_loss(v_cut, v_max, v);
      };

template <typename Param>
std::function<double(double, double, double)> retransform_loss<Param>::func
    = [](double v_cut, double v_max, double v) {
          return retransform_relative_loss(v_cut, v_max, v);
      };

namespace crosssection {
    struct ComptonKleinNishina;
}

#ifndef CROSSSECTIONDNDXINTERPOLANT_INSTANTIATION
extern template
std::function<double(double, double, double)> transform_loss<crosssection::ComptonKleinNishina>::func;
#endif

#ifndef CROSSSECTIONDNDXINTERPOLANT_INSTANTIATION
extern template
std::function<double(double, double, double)> retransform_loss<crosssection::ComptonKleinNishina>::func;
#endif

template <typename T1, typename... Args>
auto build_dndx_def(T1 const& param, ParticleDef const& p, Args... args)
{
    auto dndx = std::make_shared<CrossSectionDNDXIntegral>(param, p, args...);
    auto v_lim = AxisBuilderDNDX::v_limits { 0, 1,
        InterpolationSettings::NODES_DNDX_V };
    auto energy_lim = AxisBuilderDNDX::energy_limits();
    energy_lim.low = param.GetLowerEnergyLim(p);
    energy_lim.up = InterpolationSettings::UPPER_ENERGY_LIM;
    energy_lim.nodes = InterpolationSettings::NODES_DNDX_E ;
    auto axis_builder = AxisBuilderDNDX(v_lim, energy_lim);
    axis_builder.refine_definition_range(
        [dndx](double E) { return dndx->Calculate(E); });

    auto def = cubic_splines::BicubicSplines<double>::Definition();
    def.axis = axis_builder.Create();
    def.f = [dndx](double energy, double v) {
        auto lim = dndx->GetIntegrationLimits(energy);
        v = transform_loss<T1>::func(lim.min, lim.max, v);
        return dndx->Calculate(energy, v);
    };
    def.approx_derivates = true;
    return def;
}

class CrossSectionDNDXInterpolant : public CrossSectionDNDX {
    std::function<double(double, double, double)> transform_v;
    std::function<double(double, double, double)> retransform_v;
    cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>
        interpolant;

    std::string gen_path() const;
    std::string gen_name() const;
    size_t gen_hash(size_t) const;

public:
    template <typename Param, typename Target>
    CrossSectionDNDXInterpolant(Param param, ParticleDef const& p,
        Target const& t, std::shared_ptr<const EnergyCutSettings> cut,
        size_t hash = 0)
        : CrossSectionDNDX(param, p, t, cut, gen_hash(hash))
        , transform_v(transform_loss<Param>::func)
        , retransform_v(retransform_loss<Param>::func)
        , interpolant(build_dndx_def(param, p, t, cut), gen_path(), gen_name())
    {
        lower_energy_lim
            = interpolant.GetDefinition().GetAxis().at(0)->GetLow();
    }

    double Calculate(double E) final;

    double Calculate(double E, double v) final;

    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
