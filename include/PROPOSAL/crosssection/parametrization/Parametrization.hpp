#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace crosssection {
    struct IonizBergerSeltzerBhabha;
    template <>
    struct is_component_wise<IonizBergerSeltzerBhabha> : std::false_type {
    };

    class IonizBetheBlochRossi;
    template <>
    struct is_component_wise<IonizBetheBlochRossi> : std::false_type {
    };

    struct IonizBergerSeltzerMoller;
    template <>
    struct is_component_wise<IonizBergerSeltzerMoller> : std::false_type {
    };
} // namespace crosssection
} // namespace PROPOSAL
