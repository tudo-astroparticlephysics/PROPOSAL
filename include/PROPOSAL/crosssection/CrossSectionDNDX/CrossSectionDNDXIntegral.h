#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"

namespace PROPOSAL {
namespace crosssection {
    struct Ionization;
} // namespace crosssection
namespace detail {
    using dndx_integrand_t
        = std::function<double(double, double, double, double)>;
    using dndx_upper_lim_t
        = std::function<double(double, double, double, double)>;

    dndx_integrand_t define_dndx_integral(
        crosssection::Parametrization<Medium> const&, ParticleDef, Medium);

    dndx_integrand_t define_dndx_integral(
        crosssection::Parametrization<Component> const&, ParticleDef,
        Component);

    dndx_integrand_t define_dndx_integral(
        crosssection::Ionization const&, ParticleDef, Medium);

    dndx_upper_lim_t define_dndx_upper_lim(
        crosssection::Parametrization<Medium> const&, ParticleDef, Medium);

    dndx_upper_lim_t define_dndx_upper_lim(
        crosssection::Parametrization<Component> const&, ParticleDef,
        Component);
} // namespace detail
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDNDXIntegral : public CrossSectionDNDX {
    detail::dndx_integrand_t dndx_integral;
    detail::dndx_upper_lim_t dndx_upper_limit;

public:
    template <typename Param, typename Target>
    CrossSectionDNDXIntegral(Param param, ParticleDef const& p, Target const& t,
        std::shared_ptr<const EnergyCutSettings> cut, size_t hash = 0)
        : CrossSectionDNDX(param, p, t, cut, hash)
        , dndx_integral(detail::define_dndx_integral(param, p, t))
        , dndx_upper_limit(detail::define_dndx_upper_lim(param, p, t))
    {
    }

    double Calculate(double energy) final;

    double Calculate(double energy, double v) final;

    double GetUpperLimit(double energy, double rate) final;
};
} // namespace PROPOSAL
