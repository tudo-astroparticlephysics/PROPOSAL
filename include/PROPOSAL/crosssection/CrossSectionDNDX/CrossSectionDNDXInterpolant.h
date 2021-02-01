#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

template <typename T1, typename T2, typename... Args>
auto build_dndx_def(T1 const& param, T2 const& p_def, Args... args)
{
    auto dndx
        = std::make_shared<CrossSectionDNDXIntegral>(param, p_def, args...);
    auto def = cubic_splines::BicubicSplines<double>::Definition();
    def.axis[0] = std::make_unique<cubic_splines::ExpAxis<double>>(
        param.GetLowerEnergyLim(p_def), 1.e14, (size_t)100);
    def.axis[1]
        = std::make_unique<cubic_splines::LinAxis<double>>(0., 1., (size_t)100);
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

    std::string gen_name()
    {
        return std::string("dndx_") + std::to_string(GetHash())
            + std::string(".txt");
    }

public:
    template <typename... Args>
    CrossSectionDNDXInterpolant(Args... args)
        : CrossSectionDNDX(args...)
        , interpolant(build_dndx_def(args...), "/tmp", gen_name())
    {
        /* auto logger = Logging::Get("PROPOSAL.CrossSectionDEDX"); */
        /* logger->debug("Interpolationtables successfully build."); */
    }

    double Calculate(double E) final;

    double Calculate(double E, double v) final;

    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
