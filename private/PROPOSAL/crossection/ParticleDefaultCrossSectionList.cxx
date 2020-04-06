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

#include <unordered_map>

namespace PROPOSAL {

InterpolationDef std_interpolation_def;

CrossSectionList BuildEMinusStdCrossSections(std::shared_ptr<Medium> medium,
    std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& electron ) noexcept
{
    BremsKelnerKokoulinPetrukhin brems(electron, medium, 1.0, true);
    EpairKelnerKokoulinPetrukhin epair(electron, medium, 1.0, true);
    IonizBetheBlochRossi ioniz(electron, medium, cut, 1.0);
    ShadowButkevichMikhailov shadow;
    PhotoAbramowiczLevinLevyMaor97 photo(electron, medium, 1.0, shadow);

    CrossSectionList crosss;
    crosss.emplace_back(
        std::make_shared<BremsInterpolant>(brems, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<EpairInterpolant>(epair, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<IonizInterpolant>(ioniz, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<PhotoInterpolant>(photo, cut, std_interpolation_def));

    return crosss;
};

CrossSectionList BuildMuMinusStdCrossSections(std::shared_ptr<Medium> medium,
    std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& muon) noexcept
{
    BremsKelnerKokoulinPetrukhin brems(muon, medium, 1.0, true);
    EpairKelnerKokoulinPetrukhin epair(muon, medium, 1.0, true);
    IonizBetheBlochRossi ioniz(muon, medium, cut, 1.0);
    ShadowButkevichMikhailov shadow;
    PhotoAbramowiczLevinLevyMaor97 photo(muon, medium, 1.0, shadow);

    CrossSectionList crosss;

    crosss.emplace_back(
        std::make_shared<BremsInterpolant>(brems, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<EpairInterpolant>(epair, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<IonizInterpolant>(ioniz, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<PhotoInterpolant>(photo, cut, std_interpolation_def));

    return crosss;
};

CrossSectionList BuildTauMinusStdCrossSections(std::shared_ptr<Medium> medium,
    std::shared_ptr<EnergyCutSettings> cut,
    const ParticleDef& tau) noexcept
{
    BremsKelnerKokoulinPetrukhin brems(tau, medium, 1.0, true);
    EpairKelnerKokoulinPetrukhin epair(tau, medium, 1.0, true);
    IonizBetheBlochRossi ioniz(tau, medium, cut, 1.0);
    ShadowButkevichMikhailov shadow;
    PhotoAbramowiczLevinLevyMaor97 photo(tau, medium, 1.0, shadow);

    CrossSectionList crosss;

    crosss.emplace_back(
        std::make_shared<BremsInterpolant>(brems, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<EpairInterpolant>(epair, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<IonizInterpolant>(ioniz, cut, std_interpolation_def));
    crosss.emplace_back(
        std::make_shared<PhotoInterpolant>(photo, cut, std_interpolation_def));

    return crosss;
};

using CrossSectionListBuilder = CrossSectionList (*)(std::shared_ptr<Medium>,
    std::shared_ptr<EnergyCutSettings>,
    const ParticleDef&);

std::unordered_map<ParticleType, CrossSectionListBuilder>
    BuilderStdCrossSections{
        { ParticleType::EMinus, BuildEMinusStdCrossSections },
        { ParticleType::MuMinus, BuildMuMinusStdCrossSections },
        { ParticleType::TauMinus, BuildTauMinusStdCrossSections } };

CrossSectionList GetStdCrossSections(std::shared_ptr<Medium> medium,
    std::shared_ptr<EnergyCutSettings> cut,
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
