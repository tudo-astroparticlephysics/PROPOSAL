#pragma once
#include <type_traits>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crosssection/CrossSectionIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionInterpolant.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"

namespace PROPOSAL {

extern InterpolationDef std_interpolation_def;

template<typename CrossVec>
CrossVec append_cross(CrossVec& cross_vec) {
    return cross_vec;
}

enum PARAMETRIZATION { PARAM, PARTICLE, MEDIUM, CUT, INTERPOLATE };
template<typename CrossVec, typename P, typename... Args>
void append_cross(CrossVec& cross_vec, P param, Args... args) {
    cross_vec.push_back(make_crosssection(
                    get<PARAM>(param), get<PARTICLE>(param), get<MEDIUM>(param),
                    get<CUT>(param), get<INTERPOLATE>(param)
        )
    );
    append_cross(cross_vec, args...);
}

template <typename ParticleType>
struct DefaultCrossSections{
    template <typename CrossVec, typename P, typename M>
    static void Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate);

    template <typename M, typename... Args>
    static auto Get(ParticleType const& particle, M const& medium, Args... args)
    {
        auto cross = crosssection_list_t<ParticleType, M>();
        DefaultCrossSections<ParticleType>::Append(cross, particle, medium, args...);
        return cross;
    }

    template <typename TypeCheck, typename... Args, typename = std::enable_if_t<std::is_same<TypeCheck, std::false_type>::value>>
    static auto Get(ParticleType const& particle, Medium const& medium, Args... args)
    {
        auto cross = crosssection_list_t<ParticleDef, Medium>();
        DefaultCrossSections<ParticleType>::Append(cross, reinterpret_cast<ParticleDef const&>(particle), medium, args...);
        return cross;
    }
};

template <typename P, typename M, typename... Args>
auto GetStdCrossSections(P const& particle, M const& medium, Args... args)
{
    return DefaultCrossSections<P>::Get(particle, medium, args...);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<GammaDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto photopair = make_tuple(crosssection::PhotoPairTsai{}, p, m, nullptr, interpolate);
    auto compton = make_tuple(crosssection::ComptonKleinNishina{}, p, m, cut, interpolate);
    append_cross(cross_vec, photopair, compton);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<EMinusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = make_tuple(crosssection::BremsElectronScreening{ false }, p, m, cut, interpolate);
    auto epair = make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikhailov>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<EPlusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = make_tuple(crosssection::BremsElectronScreening{ false }, p, m, cut, interpolate);
    auto epair = make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikhailov>() }, p, m, cut, interpolate);
    auto annih = make_tuple(crosssection::AnnihilationHeitler{}, p, m, nullptr, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo, annih);
}


template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<MuMinusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikhailov>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<MuPlusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikhailov>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<TauMinusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikhailov>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<TauPlusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikhailov>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

} // namespace PROPOSAL
