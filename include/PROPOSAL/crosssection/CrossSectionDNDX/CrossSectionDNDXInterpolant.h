#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

class AxisBuilder_dNdx {
    using axis_t = cubic_splines::ExpAxis<double>;
    double low, up;
    size_t n;

public:
    template <typename Param, typename... Args>
    AxisBuilder_dNdx(Param const& _param, ParticleDef const& _p, Args...)
        : low(_param.GetLowerEnergyLim(_p))
        , up(UPPER_ENERGY_LIM_DEFAULT)
        , n(NODES_DNDX_E_DEFAULT)
    {
        if (UPPER_ENERGY_LIM)
            up = *UPPER_ENERGY_LIM;
        if (NODES_DNDX_E)
            n = *NODES_DNDX_E;
    }

    template <typename T1> auto refine_definition_range(T1 func)
    {
        auto i = 0u;
        auto ax = axis_t(low, up, n);
        while (not(func(ax.back_transform(i)) > 0) and i < n)
            ++i;
        if (i==n)
            throw std::logic_error("No positive values to build dNdx tables!");
        low = ax.back_transform(i);
    }

    auto Create() const { return std::make_unique<axis_t>(low, up, n); }
};

template <typename T1, typename... Args>
auto build_dndx_def(T1 const& param, ParticleDef const& p_def, Args... args)
{
    auto dndx
        = std::make_shared<CrossSectionDNDXIntegral>(param, p_def, args...);
    auto def = cubic_splines::BicubicSplines<double>::Definition();

    auto ax_energy = AxisBuilder_dNdx(param, p_def, args...);
    ax_energy.refine_definition_range(
        [dndx](double E) { return dndx->Calculate(E); });
    def.axis[0] = ax_energy.Create();

    auto nodes_v = NODES_DNDX_V_DEFAULT;
    if (NODES_DNDX_V)
        nodes_v = *NODES_DNDX_V;
    def.axis[1] = std::make_unique<cubic_splines::LinAxis<double>>(
        0., 0.999, (size_t)nodes_v);

    def.f = [dndx](double energy, double v) {
        auto lim = dndx->GetIntegrationLimits(energy);
        v = transform_relativ_loss(
            lim[CrossSectionDNDX::MIN], lim[CrossSectionDNDX::MAX], v);
        return dndx->Calculate(energy, v);
    };
    def.approx_derivates = true;
    return def;
}

class CrossSectionDNDXInterpolant : public CrossSectionDNDX {
    cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>
        interpolant;

    double lower_energy_lim;

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
        , lower_energy_lim(
              interpolant.GetDefinition().GetAxis().at(0)->GetLow())
    {
    }

    template <typename Param, typename Target>
    CrossSectionDNDXInterpolant(
        Param param, ParticleDef const& p, Target const& t, size_t hash = 0)
        : CrossSectionDNDXInterpolant(param, p, t, nullptr, hash)
    {
    }

    double Calculate(double E) final;

    double Calculate(double E, double v) final;

    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
