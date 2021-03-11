#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
inline double transform_loss_linear(double v_cut, double v_max, double v)
{
    return (v_max - v_cut) * v + v_cut;
}
double retransform_relativ_loss(double v_cut, double v_max, double v);

template <typename Param> struct transform_loss {
    static std::function<double(double, double, double)> func;
};

template <typename Param>
std::function<double(double, double, double)> transform_loss<Param>::func
    = [](double v_cut, double v_max, double v) {
          return transform_relativ_loss(v_cut, v_max, v);
      };

class ComptonKleinNishina;
template <>
std::function<double(double, double, double)> transform_loss<ComptonKleinNishina>::func;
    /* = [](double v_cut, double v_max, double v) { */
    /*       return transform_loss_linear(v_cut, v_max, v); */
    /*   }; */

template <typename T1, typename... Args>
auto build_dndx_def(T1 const& param, ParticleDef const& p, Args... args)
{
    auto dndx = std::make_shared<CrossSectionDNDXIntegral>(param, p, args...);
    auto v_lim = AxisBuilderDNDX::v_limits { 0, 1,
        InterpolationSettings::NODES_DNDX_V };
    auto energy_lim
        = AxisBuilderDNDX::energy_limits { param.GetLowerEnergyLim(p),
              InterpolationSettings::UPPER_ENERGY_LIM,
              InterpolationSettings::NODES_DNDX_E };
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
    cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>
        interpolant;

    std::string gen_path() const;
    std::string gen_name() const;

public:
    template <typename Param, typename Target>
    CrossSectionDNDXInterpolant(Param param, ParticleDef const& p,
        Target const& t, std::shared_ptr<const EnergyCutSettings> cut,
        size_t hash = 0)
        : CrossSectionDNDX(param, p, t, cut, hash)
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
