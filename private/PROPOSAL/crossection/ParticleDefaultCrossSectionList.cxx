#include "PROPOSAL/crossection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

#include "PROPOSAL/particle/ParticleDef.h"
#include <unordered_map>

using std::shared_ptr;
using std::vector;

namespace PROPOSAL {

InterpolationDef std_interpolation_def;

vector<shared_ptr<CrossSection>> BuildEMinusStdCrossSections(
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut) noexcept
{
    EMinusDef p_def;

    auto shadow = make_unique<ShadowButkevichMikhailov>();

    auto brems = make_unique<BremsKelnerKokoulinPetrukhin>(
        p_def, medium.GetComponents(), true);
    auto epair = make_unique<EpairKelnerKokoulinPetrukhin>(
        p_def, medium.GetComponents(), true);
    auto ioniz = make_unique<IonizBetheBlochRossi>(p_def, medium, *cut);
    auto photo = make_unique<PhotoAbramowiczLevinLevyMaor97>(
        p_def, medium.GetComponents(), std::move(shadow));

    auto brems_inter = std::make_shared<BremsInterpolant>(
        std::move(brems), cut, std_interpolation_def);
    auto epair_inter = std::make_shared<EpairInterpolant>(
        std::move(epair), cut, std_interpolation_def);
    auto ioniz_inter = std::make_shared<IonizInterpolant>(
        std::move(ioniz), cut, std_interpolation_def);
    auto photo_inter = std::make_shared<PhotoInterpolant>(
        std::move(photo), cut, std_interpolation_def);

    return vector<shared_ptr<CrossSection>>(
        { brems_inter, epair_inter, ioniz_inter, photo_inter });
}

vector<shared_ptr<CrossSection>> BuildMuMinusStdCrossSections(
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut) noexcept
{
    MuMinusDef p_def;

    auto shadow = make_unique<ShadowButkevichMikhailov>();

    auto brems = make_unique<BremsKelnerKokoulinPetrukhin>(
        p_def, medium.GetComponents(), true);
    auto epair = make_unique<EpairKelnerKokoulinPetrukhin>(
        p_def, medium.GetComponents(), true);
    auto ioniz = make_unique<IonizBetheBlochRossi>(p_def, medium, *cut);
    auto photo = make_unique<PhotoAbramowiczLevinLevyMaor97>(
        p_def, medium.GetComponents(), std::move(shadow));

    auto brems_inter = std::make_shared<BremsInterpolant>(
        std::move(brems), cut, std_interpolation_def);
    auto epair_inter = std::make_shared<EpairInterpolant>(
        std::move(epair), cut, std_interpolation_def);
    auto ioniz_inter = std::make_shared<IonizInterpolant>(
        std::move(ioniz), cut, std_interpolation_def);
    auto photo_inter = std::make_shared<PhotoInterpolant>(
        std::move(photo), cut, std_interpolation_def);

    return vector<shared_ptr<CrossSection>>(
        { brems_inter, epair_inter, ioniz_inter, photo_inter });
}

vector<shared_ptr<CrossSection>> BuildTauMinusStdCrossSections(
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut) noexcept
{
    TauMinusDef p_def;
    auto shadow = make_unique<ShadowButkevichMikhailov>();

    auto brems = make_unique<BremsKelnerKokoulinPetrukhin>(
        p_def, medium.GetComponents(), true);
    auto epair = make_unique<EpairKelnerKokoulinPetrukhin>(
        p_def, medium.GetComponents(), true);
    auto ioniz = make_unique<IonizBetheBlochRossi>(p_def, medium, *cut);
    auto photo = make_unique<PhotoAbramowiczLevinLevyMaor97>(
        p_def, medium.GetComponents(), std::move(shadow));

    auto brems_inter = std::make_shared<BremsInterpolant>(
        std::move(brems), cut, std_interpolation_def);
    auto epair_inter = std::make_shared<EpairInterpolant>(
        std::move(epair), cut, std_interpolation_def);
    auto ioniz_inter = std::make_shared<IonizInterpolant>(
        std::move(ioniz), cut, std_interpolation_def);
    auto photo_inter = std::make_shared<PhotoInterpolant>(
        std::move(photo), cut, std_interpolation_def);

    return vector<shared_ptr<CrossSection>>(
        { brems_inter, epair_inter, ioniz_inter, photo_inter });
}

using CrossSectionListBuilder = vector<shared_ptr<CrossSection>> (*)(
    const Medium&, std::shared_ptr<EnergyCutSettings>);

std::unordered_map<ParticleType, CrossSectionListBuilder>
    BuilderStdCrossSections{ { ParticleType::EMinus,
                                 BuildEMinusStdCrossSections },
        { ParticleType::MuMinus, BuildMuMinusStdCrossSections },
        { ParticleType::TauMinus, BuildTauMinusStdCrossSections } };

vector<shared_ptr<CrossSection>> GetStdCrossSections(const Medium& medium,
    std::shared_ptr<EnergyCutSettings> cut, const ParticleDef& p_def)
{
    auto search = BuilderStdCrossSections.find(
        static_cast<ParticleType>(p_def.particle_type));
    if (search != BuilderStdCrossSections.end())
        return (*search->second)(medium, cut);

    throw std::invalid_argument(
        "Partilce type has no standard cross section builder");
}
}
