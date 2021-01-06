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
auto build_definition(T const& param, Args... args)
{
    auto dndx = std::make_shared<CrossSectionDNDXIntegral>(param, args...);
    auto def = cubic_splines::BicubicSplines::Definition();
    def.axis[0] = std::make_unique<cubic_splines::ExpAxis>(dndx->GetLowerEnergyLim(), 1e14, (size_t)100);
    def.axis[1] = std::make_unique<cubic_splines::LinAxis>(0, 1, (size_t)100);
    def.f = [dndx](double energy, double v) {
        auto lim = dndx->GetIntegrationLimits(energy);
        v = transform_relativ_loss(std::get<CrossSectionDNDX::MIN>(lim),
            std::get<CrossSectionDNDX::MAX>(lim), v);
        return dndx->Calculate(energy, v);
    };
    return def;
}

class CrossSectionDNDXInterpolant : public CrossSectionDNDX,
                                    public CrossSectionInterpolantBase {

    cubic_splines::Interpolant<cubic_splines::BicubicSplines> interpolant;

public:
    CrossSectionDNDXInterpolant() = default;

    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDXInterpolant(Param&& _param, Particle& _p_def,
        Target& _target, std::shared_ptr<const EnergyCutSettings> _cut)
        : interpolant(build_definition(_param, _p_def, _target, _cut), "/tmp", "")
    {
    }

    double Calculate(double) final;
    double Calculate(double, double) final;
    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
