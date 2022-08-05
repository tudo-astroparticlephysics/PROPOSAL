#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace crosssection {
    struct IonizBergerSeltzerBhabha;
    template <>
    struct is_component_wise<IonizBergerSeltzerBhabha> : std::false_type {
    };

    struct IonizBetheBlochRossi;
    template <>
    struct is_component_wise<IonizBetheBlochRossi> : std::false_type {
    };

    struct IonizBergerSeltzerMoller;
    template <>
    struct is_component_wise<IonizBergerSeltzerMoller> : std::false_type {
    };

    struct AnnihilationHeitler;
    template <>
    struct is_only_stochastic<AnnihilationHeitler> : std::true_type {};

    struct PhotoeffectSauter;
    template <>
    struct is_only_stochastic<PhotoeffectSauter> : std::true_type {};

    struct PhotoMuPairBurkhardtKelnerKokoulin;
    template <>
    struct is_only_stochastic<PhotoMuPairBurkhardtKelnerKokoulin> : std::true_type {};

    struct PhotoPairTsai;
    template <>
    struct is_only_stochastic<PhotoPairTsai> : std::true_type {};

    struct PhotoPairKochMotz;
    template <>
    struct is_only_stochastic<PhotoPairKochMotz> : std::true_type {};

    struct WeakCooperSarkarMertsch;
    template <>
    struct is_only_stochastic<WeakCooperSarkarMertsch> : std::true_type {};

} // namespace crosssection
} // namespace PROPOSAL
