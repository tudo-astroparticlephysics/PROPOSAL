#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionInterpolantBase.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {
namespace crosssection {
    struct Parametrization;
}
struct ParticleDef;

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

template <typename T, typename... Args>
auto build_dndx_def(T const& param, Args... args)
{
    auto dndx = std::make_shared<CrossSectionDNDXIntegral>(param, args...);
    auto def = cubic_splines::BicubicSplines::Definition();
    def.axis[0] = std::make_unique<cubic_splines::ExpAxis>(
        dndx->GetLowerEnergyLim(), 1e14, (size_t)100);
    def.axis[1] = std::make_unique<cubic_splines::LinAxis>(0, 1, (size_t)100);
    def.f = [dndx](double energy, double v) {
        auto lim = dndx->GetIntegrationLimits(energy);
        v = transform_relativ_loss(lim[CrossSectionDNDX::MIN], lim[CrossSectionDNDX::MAX], v);
        return dndx->Calculate(energy, v);
    };
    return def;
}

class CrossSectionDNDXInterpolant : public CrossSectionDNDX {

    cubic_splines::Interpolant<cubic_splines::BicubicSplines> interpolant;

public:
    template <typename... Args>
    CrossSectionDNDXInterpolant(Args... args)
        : CrossSectionDNDX(args...)
        , interpolant(build_dndx_def(args...), "/tmp", "interpolant.txt")
    {
    }

    double Calculate(double) final;
    double Calculate(double, double) final;
    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
