#include "PROPOSAL/crossection/ParticleDefaultCrossSectionList.h"

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
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& electron) noexcept
{
    /* BremsKelnerKokoulinPetrukhin brems(electron, medium.GetComponents(),  true); */
    /* EpairKelnerKokoulinPetrukhin epair(electron, medium.GetComponents(), true); */
    /* IonizBetheBlochRossi ioniz(electron, medium, *cut); */
    /* ShadowButkevichMikhailov shadow; */
    /* PhotoAbramowiczLevinLevyMaor97 photo(electron, medium.GetComponents(), shadow); */

    vector<shared_ptr<CrossSection>> crosss;
    /* crosss.emplace_back( */
    /*     std::make_shared<BremsInterpolant>(brems, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<EpairInterpolant>(epair, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<IonizInterpolant>(ioniz, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<PhotoInterpolant>(photo, cut, std_interpolation_def)); */

    return crosss;
}

vector<shared_ptr<CrossSection>> BuildMuMinusStdCrossSections(
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& muon) noexcept
{
    /* BremsKelnerKokoulinPetrukhin brems(muon, medium.GetComponents(), true); */
    /* EpairKelnerKokoulinPetrukhin epair(muon, medium.GetComponents(), true); */
    /* IonizBetheBlochRossi ioniz(muon, medium, *cut); */
    /* ShadowButkevichMikhailov shadow; */
    /* PhotoAbramowiczLevinLevyMaor97 photo(muon, medium.GetComponents(), shadow); */

    vector<shared_ptr<CrossSection>> crosss;

    /* crosss.emplace_back( */
    /*     std::make_shared<BremsInterpolant>(brems, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<EpairInterpolant>(epair, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<IonizInterpolant>(ioniz, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<PhotoInterpolant>(photo, cut, std_interpolation_def)); */

    return crosss;
}

vector<shared_ptr<CrossSection>> BuildTauMinusStdCrossSections(
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& tau) noexcept
{
    /* BremsKelnerKokoulinPetrukhin brems(tau, medium.GetComponents(), true); */
    /* EpairKelnerKokoulinPetrukhin epair(tau, medium.GetComponents(), true); */
    /* IonizBetheBlochRossi ioniz(tau, medium, *cut); */
    /* ShadowButkevichMikhailov shadow; */
    /* PhotoAbramowiczLevinLevyMaor97 photo(tau, medium.GetComponents(), shadow); */

    vector<shared_ptr<CrossSection>> crosss;

    /* crosss.emplace_back( */
    /*     std::make_shared<BremsInterpolant>(brems, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<EpairInterpolant>(epair, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<IonizInterpolant>(ioniz, cut, std_interpolation_def)); */
    /* crosss.emplace_back( */
    /*     std::make_shared<PhotoInterpolant>(photo, cut, std_interpolation_def)); */

    return crosss;
}

using CrossSectionListBuilder
    = vector<shared_ptr<CrossSection>> (*)(const Medium&,
        std::shared_ptr<EnergyCutSettings>, const ParticleDef&);

std::unordered_map<ParticleType, CrossSectionListBuilder>
    BuilderStdCrossSections{ { ParticleType::EMinus,
                                 BuildEMinusStdCrossSections },
        { ParticleType::MuMinus, BuildMuMinusStdCrossSections },
        { ParticleType::TauMinus, BuildTauMinusStdCrossSections } };

vector<shared_ptr<CrossSection>> GetStdCrossSections(
    const Medium& medium, std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& p_def)
{
    auto search = BuilderStdCrossSections.find(
        static_cast<ParticleType>(p_def.particle_type));
    if (search != BuilderStdCrossSections.end())
        return (*search->second)(medium, cut, p_def);

    throw std::invalid_argument(
        "Partilce type has no standard cross section builder");
}
}
