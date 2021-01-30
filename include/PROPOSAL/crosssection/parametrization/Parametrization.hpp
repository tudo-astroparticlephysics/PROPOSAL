#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace crosssection {
    class IonizBergerSeltzerBhabha;
    template <>
    struct is_component_wise<IonizBergerSeltzerBhabha> : std::false_type {
    };

    class IonizBetheBlochRossi;
    template <>
    struct is_component_wise<IonizBetheBlochRossi> : std::false_type {
    };

    class IonizBergerSeltzerMoller;
    template <>
    struct is_component_wise<IonizBergerSeltzerMoller> : std::false_type {
    };
} // namespace crosssection
} // namespace PROPOSAL
