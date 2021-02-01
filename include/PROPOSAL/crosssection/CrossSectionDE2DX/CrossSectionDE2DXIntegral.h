#pragma once

#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DX.h"

#include <functional>

namespace PROPOSAL {
class Integral;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace detail {
    std::function<double(Integral&, double)> define_de2dx_integral(
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p_def, Medium const& target,
        EnergyCutSettings const& cut);

    std::function<double(Integral&, double)> define_de2dx_integral(
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p_def, Component const& target,
        EnergyCutSettings const& cut);
} // namespace detail
}

namespace PROPOSAL {
class CrossSectionDE2DXIntegral : public CrossSectionDE2DX {
    std::function<double(Integral&, double)> de2dx_integral;

public:
    template <typename... Args>
    CrossSectionDE2DXIntegral(Args... args)
        : CrossSectionDE2DX(args...)
        , de2dx_integral(detail::define_de2dx_integral(args...))
    {
    }

    double Calculate(double E) const final;
};
} // namespace PROPOSAL
