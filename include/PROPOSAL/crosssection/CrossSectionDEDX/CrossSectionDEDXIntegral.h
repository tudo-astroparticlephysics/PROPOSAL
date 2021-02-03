#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include <functional>
#include <memory>

namespace PROPOSAL {
namespace crosssection {
    struct IonizBergerSeltzerBhabha;
    struct IonizBergerSeltzerMoller;
} // namespace crosssection
class Integral;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace detail {
    using dedx_integral_t = std::function<double(Integral&, double)>;

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Component> const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Medium> const&, ParticleDef const&,
        Medium const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&);
}

class CrossSectionDEDXIntegral : public CrossSectionDEDX {
    std::function<double(Integral&, double)> dedx_integral;

public:
    template <typename... Args>
    CrossSectionDEDXIntegral(Args... args)
        : CrossSectionDEDX(args...)
        , dedx_integral(detail::define_dedx_integral(args...))
    {
    }

    double Calculate(double E) const final;
};
} // namespace PROPOSAL
