#include "PROPOSAL/crossection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

#include "PROPOSAL/particle/ParticleDef.h"
#include <unordered_map>

using std::make_shared;
using std::shared_ptr;
using std::unique_ptr;
using std::vector;

namespace PROPOSAL {
template <typename Param>
shared_ptr<CrossSection> make_crosssection(
    Param&& param, shared_ptr<const EnergyCutSettings> cuts, bool interpolate)
{
    if (interpolate)
        return make_shared<CrossSectionInterpolant<Param>>(
            param, cuts, InterpolationDef());
    return make_shared<CrossSectionIntegral<Param>>(param, cuts);
}

vector<shared_ptr<CrossSection>> BuildEMinusStdCrossSections(
    std::shared_ptr<EnergyCutSettings> cut) noexcept
{
    BremsKelnerKokoulinPetrukhin brems{ false };
    EpairKelnerKokoulinPetrukhin epair{ false };
    IonizBetheBlochRossi ioniz{ EnergyCutSettings(*cut) };
    PhotoAbramowiczLevinLevyMaor97 photo{
        make_unique<ShadowButkevichMikhailov>()
    };

    vector<shared_ptr<CrossSection>> cross_list;
    cross_list.emplace_back(make_crosssection(brems, cut, false));
    cross_list.emplace_back(make_crosssection(epair, cut, false));
    cross_list.emplace_back(make_crosssection(ioniz, cut, false));
    cross_list.emplace_back(make_crosssection(photo, cut, false));

    return cross_list;
}

/* vector<shared_ptr<CrossSection>> BuildMuMinusStdCrossSections( */
/*     std::shared_ptr<EnergyCutSettings> cut) noexcept */
/* { */
/*     MuMinusDef p_def; */

/*     auto shadow = make_unique<ShadowButkevichMikhailov>(); */

/*     auto brems = make_unique<BremsKelnerKokoulinPetrukhin>(false); */
/*     auto epair = make_unique<EpairKelnerKokoulinPetrukhin>(false); */
/*     auto ioniz = make_unique<IonizBetheBlochRossi>(*cut); */
/*     auto photo =
 * make_unique<PhotoAbramowiczLevinLevyMaor97>(std::move(shadow)); */

/*     auto brems_inter = std::make_shared<BremsInterpolant>( */
/*         std::move(brems), cut, std_interpolation_def); */
/*     auto epair_inter = std::make_shared<EpairInterpolant>( */
/*         std::move(epair), cut, std_interpolation_def); */
/*     auto ioniz_inter = std::make_shared<IonizInterpolant>( */
/*         std::move(ioniz), cut, std_interpolation_def); */
/*     auto photo_inter = std::make_shared<PhotoInterpolant>( */
/*         std::move(photo), cut, std_interpolation_def); */

/*     return vector<shared_ptr<CrossSection>>( */
/*         { brems_inter, epair_inter, ioniz_inter, photo_inter }); */
/* } */

/* vector<shared_ptr<CrossSection>> BuildTauMinusStdCrossSections( */
/*     std::shared_ptr<EnergyCutSettings> cut) noexcept */
/* { */
/*     TauMinusDef p_def; */
/*     auto shadow = make_unique<ShadowButkevichMikhailov>(); */

/*     auto brems = make_unique<BremsKelnerKokoulinPetrukhin>(false); */
/*     auto epair = make_unique<EpairKelnerKokoulinPetrukhin>(false); */
/*     auto ioniz = make_unique<IonizBetheBlochRossi>(*cut); */
/*     auto photo =
 * make_unique<PhotoAbramowiczLevinLevyMaor97>(std::move(shadow)); */

/*     auto brems_inter = std::make_shared<BremsInterpolant>( */
/*         std::move(brems), cut, std_interpolation_def); */
/*     auto epair_inter = std::make_shared<EpairInterpolant>( */
/*         std::move(epair), cut, std_interpolation_def); */
/*     auto ioniz_inter = std::make_shared<IonizInterpolant>( */
/*         std::move(ioniz), cut, std_interpolation_def); */
/*     auto photo_inter = std::make_shared<PhotoInterpolant>( */
/*         std::move(photo), cut, std_interpolation_def); */

/*     return vector<shared_ptr<CrossSection>>( */
/*         { brems_inter, epair_inter, ioniz_inter, photo_inter }); */
/* } */

using CrossSectionListBuilder
    = vector<shared_ptr<CrossSection>> (*)(std::shared_ptr<EnergyCutSettings>);

std::unordered_map<ParticleType, CrossSectionListBuilder>
    BuilderStdCrossSections{
        { ParticleType::EMinus, BuildEMinusStdCrossSections },
        /* { ParticleType::MuMinus, BuildMuMinusStdCrossSections }, */
        /* { ParticleType::TauMinus, BuildTauMinusStdCrossSections } */
    };

vector<shared_ptr<CrossSection>> GetStdCrossSections(
    std::shared_ptr<EnergyCutSettings> cut, const ParticleDef& p_def)
{
    auto search = BuilderStdCrossSections.find(
        static_cast<ParticleType>(p_def.particle_type));
    if (search != BuilderStdCrossSections.end())
        return (*search->second)(cut);

    throw std::invalid_argument(
        "Partilce type has no standard cross section builder");
}
}
