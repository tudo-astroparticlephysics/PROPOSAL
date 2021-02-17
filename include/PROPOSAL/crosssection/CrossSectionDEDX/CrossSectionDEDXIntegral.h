#pragma once

#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include <functional>

namespace PROPOSAL {
namespace crosssection {
    struct IonizBergerSeltzerBhabha;
    struct IonizBergerSeltzerMoller;
    struct IonizBetheBlochRossi;
    struct ComptonKleinNishina;
    class EpairProduction;
    class Photonuclear;
    class MupairProduction;
} // namespace crosssection
} // namespace PROPOSAL

namespace PROPOSAL {
namespace detail {
    using dedx_integral_t = std::function<double(double)>;

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Component> const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Medium> const&, ParticleDef const&,
        Medium const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha const&, ParticleDef const&,
        Medium const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller const&, ParticleDef const&,
        Medium const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBetheBlochRossi const&, ParticleDef const&,
        Medium const& medium, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::ComptonKleinNishina const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::EpairProduction const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::Photonuclear const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);

    dedx_integral_t define_dedx_integral(
        crosssection::MupairProduction const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);
} // namespace detail
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDEDXIntegral : public CrossSectionDEDX {
    std::function<double(double)> dedx_integral;

public:
    template <typename Param, typename Target>
    CrossSectionDEDXIntegral(Param const& param, ParticleDef const& p,
        Target const& t, EnergyCutSettings const& cut, size_t hash = 0)
        : CrossSectionDEDX(param, p, t, cut, hash)
        , dedx_integral(detail::define_dedx_integral(param, p, t, cut))
    {
    }

    double Calculate(double E) const final;
};
} // namespace PROPOSAL
