#pragma once
#include <type_traits>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

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
            param, p_def, medium, cuts, InterpolationDef());
    return PROPOSAL::make_unique<CrossSectionIntegral<Param, P, M>>(
        param, p_def, medium, cuts);
}

template <typename M>
crosssection_list_t<EMinusDef, M> GetStdCrossSections(const EMinusDef& e,
    M&& medium, std::shared_ptr<const EnergyCutSettings> cut)
{
    BremsKelnerKokoulinPetrukhin brems{ false };
    EpairKelnerKokoulinPetrukhin epair{ false };
    IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    PhotoAbramowiczLevinLevyMaor97 photo{
        make_unique<ShadowButkevichMikhailov>()
    };
    crosssection_list_t<EMinusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, e, medium, cut, false));
    cross.emplace_back(make_crosssection(epair, e, medium, cut, false));
    cross.emplace_back(make_crosssection(ioniz, e, medium, cut, false));
    cross.emplace_back(make_crosssection(photo, e, medium, cut, false));
    return cross;
}

template <typename M>
crosssection_list_t<MuMinusDef, M> GetStdCrossSections(const MuMinusDef& mu,
    M&& medium, std::shared_ptr<const EnergyCutSettings> cut)
{
    BremsKelnerKokoulinPetrukhin brems{ false };
    EpairKelnerKokoulinPetrukhin epair{ false };
    IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    PhotoAbramowiczLevinLevyMaor97 photo{
        make_unique<ShadowButkevichMikhailov>()
    };
    crosssection_list_t<MuMinusDef, M> cross;
    cross.emplace_back(make_crosssection(brems, mu, medium, cut, false));
    cross.emplace_back(make_crosssection(epair, mu, medium, cut, false));
    cross.emplace_back(make_crosssection(ioniz, mu, medium, cut, false));
    cross.emplace_back(make_crosssection(photo, mu, medium, cut, false));
    return cross;
}

} // namespace PROPOSAL
