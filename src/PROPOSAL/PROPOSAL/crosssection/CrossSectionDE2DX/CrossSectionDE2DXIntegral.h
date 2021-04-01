#pragma once

#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DX.h"

#include <functional>
#include <memory>

namespace PROPOSAL {
namespace crosssection {
    template <typename T> class Parametrization;
} // namespace crosssection
struct ParticleDef;
class Medium;
class Component;
class EnergyCutSettings;
class Integral;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace detail {
    std::function<double(double)> define_de2dx_integral(
        crosssection::Parametrization<Medium> const&, ParticleDef const&,
        Medium const&, EnergyCutSettings const&);

    std::function<double(double)> define_de2dx_integral(
        crosssection::Parametrization<Component> const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);
} // namespace detail
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDE2DXIntegral : public CrossSectionDE2DX {
    std::function<double(double)> de2dx_integral;

public:
    template <typename Param, typename Target>
    CrossSectionDE2DXIntegral(Param const& param, ParticleDef const& p,
        Target const& t, EnergyCutSettings const& cut, size_t hash = 0)
        : CrossSectionDE2DX(param, p, t, cut, hash)
        , de2dx_integral(detail::define_de2dx_integral(param, p, t, cut))
    {
    }

    double Calculate(double E) const final;
};
} // namespace PROPOSAL
