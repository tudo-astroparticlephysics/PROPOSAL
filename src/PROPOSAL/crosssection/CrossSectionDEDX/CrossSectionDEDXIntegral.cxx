#include "PROPOSAL/crosssection/parametrization/Ionization.h"

namespace PROPOSAL {
namespace detail {
    template <>
    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&)
    {
        return [param, p_def, medium](Integral&, double E) {
            return param.FunctionToDEdxIntegral(p_def, medium, E, 0.) / E;
        };
    }

    template <>
    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller param, ParticleDef const& p_def,
        Medium const& medium, EnergyCutSettings const&)
    {
        return [param, p_def, medium](Integral&, double E) {
            return param.FunctionToDEdxIntegral(p_def, medium, E, 0.) / E;
        };
    }
}
}
