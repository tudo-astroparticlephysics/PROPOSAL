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

using std::decay;
using std::make_shared;
using std::unique_ptr;
using std::vector;

namespace PROPOSAL {

extern InterpolationDef std_interpolation_def;

template <typename Param, typename P, typename M>
unique_ptr<crosssection_t<P, M>> make_crosssection(Param&& param, P&& p_def,
    M&& medium, shared_ptr<const EnergyCutSettings> cuts, bool interpolate)
{
    if (interpolate)
        return PROPOSAL::make_unique<CrossSectionInterpolant<Param, P, M>>(
            param, std::forward<P>(p_def), std::forward<M>(medium), cuts, InterpolationDef());
    return PROPOSAL::make_unique<CrossSectionIntegral<Param, P, M>>(
        param, std::forward<P>(p_def), std::forward<M>(medium), cuts);
}
template <typename P, typename M>
crosssection_list_t<EMinusDef, M> GetStdCrossSections(P particle,
    M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    std::cout << "name: " << particle.name << std::endl;
    throw std::logic_error("No default crosssection for this particle provided.");
}

template <typename M>
crosssection_list_t<GammaDef, M> GetStdCrossSections(const GammaDef& gamma,
                                                      M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::PhotoPairTsai photopair {};
    crosssection::ComptonKleinNishina compton {};
    crosssection_list_t<GammaDef, M> cross;
    cross.emplace_back(make_crosssection(photopair, gamma, medium, nullptr, interpolate));
    cross.emplace_back(make_crosssection(compton, gamma, medium, cut, interpolate));

    return cross;
}

template <typename M>
crosssection_list_t<EMinusDef, M> GetStdCrossSections(const EMinusDef& e,
    M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::BremsElectronScreening brems{ false };
    crosssection::EpairKelnerKokoulinPetrukhin epair{ false };
    crosssection::IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    crosssection::PhotoAbramowiczLevinLevyMaor97 photo{
        make_unique<crosssection::ShadowButkevichMikhailov>()
    };
    crosssection_list_t<EMinusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(epair, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(ioniz, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(photo, e, medium, cut, interpolate));
    return cross;
}

template <typename M>
crosssection_list_t<EPlusDef, M> GetStdCrossSections(const EPlusDef& e,
    M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::BremsElectronScreening brems{ false };
    crosssection::EpairKelnerKokoulinPetrukhin epair{ false };
    crosssection::IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    crosssection::PhotoAbramowiczLevinLevyMaor97 photo{
        make_unique<crosssection::ShadowButkevichMikhailov>()
    };
    crosssection::AnnihilationHeitler annih {};
    crosssection_list_t<EPlusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(epair, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(ioniz, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(photo, e, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(annih, e, medium, nullptr, interpolate));
    return cross;
}

template <typename M>
crosssection_list_t<MuMinusDef, M> GetStdCrossSections(const MuMinusDef& mu,
    M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::BremsKelnerKokoulinPetrukhin brems{ false };
    crosssection::EpairKelnerKokoulinPetrukhin epair{ false };
    crosssection::IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    crosssection::PhotoAbramowiczLevinLevyMaor97 photo{
        make_unique<crosssection::ShadowButkevichMikhailov>()
    };
    crosssection_list_t<MuMinusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, mu, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(epair, mu, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(ioniz, mu, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(photo, mu, medium, cut, interpolate));
    return cross;
}

template <typename M>
crosssection_list_t<MuPlusDef, M> GetStdCrossSections(const MuPlusDef& mu,
                                                       M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::BremsKelnerKokoulinPetrukhin brems{ false };
    crosssection::EpairKelnerKokoulinPetrukhin epair{ false };
    crosssection::IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    crosssection::PhotoAbramowiczLevinLevyMaor97 photo{
            make_unique<crosssection::ShadowButkevichMikhailov>()
    };
    crosssection_list_t<MuPlusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, mu, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(epair, mu, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(ioniz, mu, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(photo, mu, medium, cut, interpolate));
    return cross;
}

template <typename M>
crosssection_list_t<TauMinusDef, M> GetStdCrossSections(const TauMinusDef& tau,
                                                       M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::BremsKelnerKokoulinPetrukhin brems{ false };
    crosssection::EpairKelnerKokoulinPetrukhin epair{ false };
    crosssection::IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    crosssection::PhotoAbramowiczLevinLevyMaor97 photo{
            make_unique<crosssection::ShadowButkevichMikhailov>()
    };
    crosssection_list_t<TauMinusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, tau, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(epair, tau, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(ioniz, tau, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(photo, tau, medium, cut, interpolate));
    return cross;
}

template <typename M>
crosssection_list_t<TauPlusDef, M> GetStdCrossSections(const TauPlusDef& tau,
                                                        M&& medium, std::shared_ptr<const EnergyCutSettings> cut, bool interpolate=true)
{
    crosssection::BremsKelnerKokoulinPetrukhin brems{ false };
    crosssection::EpairKelnerKokoulinPetrukhin epair{ false };
    crosssection::IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    crosssection::PhotoAbramowiczLevinLevyMaor97 photo{
            make_unique<crosssection::ShadowButkevichMikhailov>()
    };
    crosssection_list_t<TauPlusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, tau, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(epair, tau, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(ioniz, tau, medium, cut, interpolate));
    cross.emplace_back(make_crosssection(photo, tau, medium, cut, interpolate));
    return cross;
}

} // namespace PROPOSAL
