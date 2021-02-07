#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

template <typename T1, typename... Args>
auto build_dndx_def(T1 const& param, ParticleDef const& p, Args... args)
{
    auto dndx = std::make_shared<CrossSectionDNDXIntegral>(param, p, args...);
    auto v_lim = AxisBuilderDNDX::v_limits { 0, 1, NODES_DNDX_V_DEFAULT };
    auto energy_lim
        = AxisBuilderDNDX::energy_limits { param.GetLowerEnergyLim(p),
              UPPER_ENERGY_LIM_DEFAULT, NODES_DNDX_E_DEFAULT };
    auto axis_builder = AxisBuilderDNDX(v_lim, energy_lim);
    axis_builder.refine_definition_range(
        [dndx](double E) { return dndx->Calculate(E); });

    auto def = cubic_splines::BicubicSplines<double>::Definition();
    def.axis = axis_builder.Create();
    def.f = [dndx](double energy, double v) {
        auto lim = dndx->GetIntegrationLimits(energy);
        v = transform_relativ_loss(lim.min, lim.max, v);
        return dndx->Calculate(energy, v);
    };
    def.approx_derivates = true;
    return def;
}

class CrossSectionDNDXInterpolant : public CrossSectionDNDX {
    cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>
        interpolant;

    std::string gen_name()
    {
        return std::string("dndx_") + std::to_string(GetHash())
            + std::string(".txt");
    }

public:
    template <typename Param, typename Target>
    CrossSectionDNDXInterpolant(Param param, ParticleDef const& p,
        Target const& t, std::shared_ptr<const EnergyCutSettings> cut,
        size_t hash = 0)
        : CrossSectionDNDX(param, p, t, cut, hash)
        , interpolant(build_dndx_def(param, p, t, cut), "/tmp", gen_name())
    {
        lower_energy_lim
            = interpolant.GetDefinition().GetAxis().at(0)->GetLow();
    }

    double Calculate(double E) final;

    double Calculate(double E, double v) final;

    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
