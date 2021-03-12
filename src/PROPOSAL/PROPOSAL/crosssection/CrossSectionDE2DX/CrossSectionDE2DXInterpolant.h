
#pragma once

#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDE2DX/AxisBuilderDE2DX.h"

#include <type_traits>

#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

template <typename T1, typename T2, typename... Args>
auto build_de2dx_def(T1 const& param, T2 const& p_def, Args... args)
{
    auto de2dx
        = std::make_shared<CrossSectionDE2DXIntegral>(param, p_def, args...);
    auto ax = AxisBuilderDE2DX(param.GetLowerEnergyLim(p_def));
    ax.refine_definition_range([de2dx](double E) { return de2dx->Calculate(E); });

    auto def = cubic_splines::CubicSplines<double>::Definition();
    def.f = [de2dx](double E) { return de2dx->Calculate(E); };
    def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1., 0.);
    def.axis = ax.Create();
    return def;
}

class CrossSectionDE2DXInterpolant : public CrossSectionDE2DX {

    cubic_splines::Interpolant<cubic_splines::CubicSplines<double>> interpolant;

    std::string gen_name() const;
    std::string gen_path() const;
    size_t gen_hash(size_t) const;

public:
    template <typename Param, typename Target>
    CrossSectionDE2DXInterpolant(Param const& param, ParticleDef const& p,
        Target const& t, EnergyCutSettings const& cut, size_t hash = 0)
        : CrossSectionDE2DX(param, p, t, cut, gen_hash(hash))
        , interpolant(build_de2dx_def(param, p, t, cut), gen_path(), gen_name())
    {
            lower_energy_lim = interpolant.GetDefinition().GetAxis().GetLow();
    }

    double Calculate(double E) const final;
};
} // namespace PROPOSAL
