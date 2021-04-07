
#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/AxisBuilderDEDX.h"
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"

namespace PROPOSAL {

template <typename T1, typename... Args>
auto build_dedx_def(T1 const& param, ParticleDef const& p, Args... args)
{
    auto dedx = std::make_shared<CrossSectionDEDXIntegral>(param, p, args...);
    auto ax = AxisBuilderDEDX(param.GetLowerEnergyLim(p));
    ax.refine_definition_range([dedx](double E) { return dedx->Calculate(E); });

    auto def = cubic_splines::CubicSplines<double>::Definition();
    def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1., 0.);
    def.f = [dedx](double E) { return dedx->Calculate(E); };
    def.axis = ax.Create();
    return def;
}

class CrossSectionDEDXInterpolant : public CrossSectionDEDX {
    cubic_splines::Interpolant<cubic_splines::CubicSplines<double>> interpolant;

    std::string gen_path() const;
    std::string gen_name() const;
    size_t gen_hash(size_t) const;

public:
    template <typename Param, typename Target>
    CrossSectionDEDXInterpolant(Param const& param, ParticleDef const& p,
        Target const& t, EnergyCutSettings const& cut, size_t hash = 0)
        : CrossSectionDEDX(param, p, t, cut, gen_hash(hash))
        , interpolant(build_dedx_def(param, p, t, cut), gen_path(), gen_name())
    {
        lower_energy_lim = interpolant.GetDefinition().GetAxis().GetLow();
    }

    double Calculate(double E) const final;
};
} // namespace PROPOSAL
