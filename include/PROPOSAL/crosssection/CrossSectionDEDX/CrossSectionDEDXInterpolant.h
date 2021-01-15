
#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

template <typename T1, typename T2, typename... Args>
auto build_dedx_def(T1 const& param, T2 const& p_def, Args... args)
{
    auto dedx
        = std::make_shared<CrossSectionDEDXIntegral>(param, p_def, args...);
    auto def = cubic_splines::CubicSplines::Definition();
    def.axis = std::make_unique<cubic_splines::ExpAxis>(
        param.GetLowerEnergyLim(p_def), 1e14, (size_t)100);
    def.f = [dedx](double E) { return dedx->Calculate(E); };
    return def;
}

class CrossSectionDEDXInterpolant : public CrossSectionDEDX {

    cubic_splines::Interpolant<cubic_splines::CubicSplines> interpolant;

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
    {
    }

    double Calculate(double E) final;
};
} // namespace PROPOSAL
