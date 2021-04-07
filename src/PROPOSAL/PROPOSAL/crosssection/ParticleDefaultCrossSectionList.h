#pragma once
#include <type_traits>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crosssection/CrossSectionIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionInterpolant.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"

#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"


namespace PROPOSAL {

template<typename CrossVec>
CrossVec append_cross(CrossVec& cross_vec) {
    return cross_vec;
}

enum PARAMETRIZATION { PARAM, PARTICLE, MEDIUM, CUT, INTERPOLATE };
template<typename CrossVec, typename P, typename... Args>
void append_cross(CrossVec& cross_vec, P param, Args... args) {
    cross_vec.push_back(make_crosssection(
                    std::get<PARAM>(param), std::get<PARTICLE>(param), std::get<MEDIUM>(param),
                    std::get<CUT>(param), std::get<INTERPOLATE>(param)
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
        auto cross = std::vector<std::shared_ptr<CrossSectionBase>>();
        DefaultCrossSections<ParticleType>::Append(cross, particle, medium, args...);
        return cross;
    }

    template <typename M, typename... Args>
    static auto Get(ParticleDef const& particle, M const& medium, Args... args)
    {
        auto cross = std::vector<std::shared_ptr<CrossSectionBase>>();
        DefaultCrossSections<ParticleType>::Append(cross, particle, medium, args...);
        return cross;
    }
};

template <typename P, typename M, typename... Args>
auto GetStdCrossSections(P const& particle, M const& medium, Args... args)
{
    return DefaultCrossSections<P>::Get(particle, medium, args...);
}

template <typename M, typename... Args>
auto GetStdCrossSections(ParticleDef const& particle, M const& medium, Args... args)
{
    switch (particle.particle_type) {
        case 11:
            return DefaultCrossSections<EMinusDef>::Get(particle, medium, args...);
        case -11:
            return DefaultCrossSections<EPlusDef>::Get(particle, medium, args...);
        case 13:
            return DefaultCrossSections<MuMinusDef>::Get(particle, medium, args...);
        case -13:
            return DefaultCrossSections<MuPlusDef>::Get(particle, medium, args...);
        case 15:
            return DefaultCrossSections<TauMinusDef>::Get(particle, medium, args...);
        case -15:
            return DefaultCrossSections<TauPlusDef>::Get(particle, medium, args...);
        case 22:
            return DefaultCrossSections<GammaDef>::Get(particle, medium, args...);
        default:
            throw std::invalid_argument("No StdCrossSection found for particle_type " + std::to_string(particle.particle_type));
    }

    //return DefaultCrossSections<P>::Get(particle, medium, args...);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<GammaDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto photopair = std::make_tuple(crosssection::PhotoPairTsai{}, p, m, nullptr, interpolate);
    auto compton = std::make_tuple(crosssection::ComptonKleinNishina{}, p, m, cut, interpolate);
    append_cross(cross_vec, photopair, compton);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<EMinusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = std::make_tuple(crosssection::BremsElectronScreening{ false }, p, m, cut, interpolate);
    auto epair = std::make_tuple(crosssection::EpairForElectronPositron{ false }, p, m, cut, interpolate);
    auto ioniz = std::make_tuple(crosssection::IonizBergerSeltzerBhabha{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = std::make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<EPlusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = std::make_tuple(crosssection::BremsElectronScreening{ false }, p, m, cut, interpolate);
    auto epair = std::make_tuple(crosssection::EpairForElectronPositron{ false }, p, m, cut, interpolate);
    auto ioniz = std::make_tuple(crosssection::IonizBergerSeltzerMoller{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = std::make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() }, p, m, cut, interpolate);
    auto annih = std::make_tuple(crosssection::AnnihilationHeitler{}, p, m, nullptr, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo, annih);
}


template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<MuMinusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = std::make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = std::make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = std::make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = std::make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<MuPlusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = std::make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = std::make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = std::make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = std::make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<TauMinusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = std::make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = std::make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = std::make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = std::make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

template<>
template <typename CrossVec, typename P, typename M>
void DefaultCrossSections<TauPlusDef>::Append(CrossVec& cross_vec, P p, M m,std::shared_ptr<const EnergyCutSettings> cut, bool interpolate)
{
    auto brems = std::make_tuple(crosssection::BremsKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto epair = std::make_tuple(crosssection::EpairKelnerKokoulinPetrukhin{ false }, p, m, cut, interpolate);
    auto ioniz = std::make_tuple(crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cut) }, p, m, cut, interpolate);
    auto photo = std::make_tuple(crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() }, p, m, cut, interpolate);
    append_cross(cross_vec, brems, epair, ioniz, photo);
}

} // namespace PROPOSAL
