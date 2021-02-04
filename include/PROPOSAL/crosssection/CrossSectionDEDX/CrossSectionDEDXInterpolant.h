
#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <type_traits>

#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

class AxisBuilder_dEdx {
    using axis_t = cubic_splines::ExpAxis<double>;
    double low, up;
    size_t n;

public:
    template <typename Param, typename... Args>
    AxisBuilder_dEdx(Param const& _param, ParticleDef const& _p, Args...)
        : low(_param.GetLowerEnergyLim(_p))
        , up(UPPER_ENERGY_LIM_DEFAULT)
        , n(NODES_DEDX_DEFAULT)
    {
        if (UPPER_ENERGY_LIM)
            up = *UPPER_ENERGY_LIM;
        if (NODES_DEDX)
            n = *NODES_DEDX;
    }

    template <typename T1> auto refine_definition_range(T1 func)
    {
        auto i = 0u;
        auto ax = axis_t(low, up, n);
        while (not(func(ax.back_transform(i)) > 0))
            ++i;
        low = ax.back_transform(i);
    }

    auto Create() const { return std::make_unique<axis_t>(low, up, n); }
};

template <typename T1, typename... Args>
auto build_dedx_def(T1 const& param, Args... args)
{
    auto dedx = std::make_shared<CrossSectionDEDXIntegral>(param, args...);
    auto ax = AxisBuilder_dEdx(param, args...);
    ax.refine_definition_range([dedx](double E) { return dedx->Calculate(E); });

    auto def = cubic_splines::CubicSplines<double>::Definition();
    def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1., 0.);
    def.f = [dedx](double E) { return dedx->Calculate(E); };
    def.axis = ax.Create();
    return def;
}

class CrossSectionDEDXInterpolant : public CrossSectionDEDX {
    cubic_splines::Interpolant<cubic_splines::CubicSplines<double>> interpolant;

    double lower_energy_lim;

    std::string gen_name()
    {
        return std::string("dedx_") + std::to_string(GetHash())
            + std::string(".txt");
    }

public:
    template <typename... Args>
    CrossSectionDEDXInterpolant(Args... args)
        : CrossSectionDEDX(args...)
        , interpolant(build_dedx_def(args...), "/tmp", gen_name())
        , lower_energy_lim(interpolant.GetDefinition().GetAxis().GetLow())
    {
        /* logger->debug("Interpolationtables successfully build."); */
    }

    double Calculate(double E) const final {
        if (E < lower_energy_lim)
            return 0.;
        return interpolant.evaluate(E);
    }
};
} // namespace PROPOSAL
