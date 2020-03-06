
#include <PROPOSAL/crossection/factories/PhotoPairFactory.h>
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, Utility::Definition const& util_definition)
{
    std::stringstream ss;
    ss << " Utility Definition (" << &util_definition << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Annihilation Definition:\n" << util_definition.annihilation_def << std::endl;
    os << "Bremsstrahlung Definition:\n" << util_definition.brems_def << std::endl;
    os << "Compton Definition:\n" << util_definition.compton_def << std::endl;
    os << "EPair Production Definition:\n" << util_definition.epair_def << std::endl;
    os << "Ionization Definition:\n" << util_definition.ioniz_def << std::endl;
    os << "MuPair Production Definition:\n" << util_definition.mupair_def << std::endl;
    os << "Photonuclear Definition:\n" << util_definition.photo_def << std::endl;
    os << "PhotoPair Production Definition:\n" << util_definition.photopair_def << std::endl;
    os << "Weak Interaction Definition:\n" << util_definition.weak_def << std::endl;

    os << Helper::Centered(60, "");
    return os;
}

} // namespace PROPOSAL
/******************************************************************************
 *                            Propagation utility                              *
 ******************************************************************************/

Utility::Definition::Definition()
    // : do_interpolation(true)
    // , interpolation_def()
    : brems_def()
    , compton_def()
    , photo_def()
    , epair_def()
    , ioniz_def()
    , mupair_def()
    , weak_def()
    , photopair_def()
    , annihilation_def()
{
}

Utility::Definition::Definition(const nlohmann::json& config)
    // : do_interpolation(true)
    // , interpolation_def()
    : brems_def()
    , compton_def()
    , photo_def()
    , epair_def()
    , ioniz_def()
    , mupair_def()
    , weak_def()
    , photopair_def()
    , annihilation_def()
{
    assert(config.is_object());

    std::string brems_str = config.value("brems", "none");
    brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(brems_str);
    brems_def.multiplier = config.value("brems_multiplier", 1.0);
    brems_def.lpm_effect = config.value("lpm", true);

    std::string compton_str = config.value("compton", "none");
    compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(compton_str);
    compton_def.multiplier = config.value("compton_multiplier", 1.0);

    std::string photo_str = config.value("photo", "none");
    photo_def.parametrization = PhotonuclearFactory::Get().GetEnumFromString(photo_str);
    photo_def.multiplier = config.value("photo_multiplier", 1.0);
    std::string photo_shadow = config.value("photo_shadow", "shadow_none");
    photo_def.shadow = PhotonuclearFactory::Get().GetShadowEnumFromString(photo_shadow);
    photo_def.hard_component = config.value("photo_hard_component", true);

    std::string epair_str = config.value("epair", "none");
    epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(epair_str);
    epair_def.multiplier = config.value("epair_multiplier", 1.0);
    epair_def.lpm_effect = config.value("lpm", true);

    std::string ioniz_str = config.value("ioniz", "none");
    ioniz_def.parametrization = IonizationFactory::Get().GetEnumFromString(ioniz_str);
    ioniz_def.multiplier = config.value("ioniz_multiplier", 1.0);

    std::string mupair_str = config.value("mupair", "none");
    mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(mupair_str);
    mupair_def.multiplier = config.value("mupair_multiplier", 1.0);
    mupair_def.particle_output = config.value("mupair_particle_output", true);

    std::string weak_str = config.value("weak", "none");
    weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(weak_str);
    weak_def.multiplier = config.value("weak_multiplier", 1.0);

    std::string photopair_str = config.value("photopair", "none");
    photopair_def.parametrization = PhotoPairFactory::Get().GetEnumFromString(photopair_str);
    photopair_def.multiplier = config.value("photopair_multiplier", 1.0);
    std::string photoangle_str = config.value("photangle", "PhotoAngleNoDeflection");
    photopair_def.photoangle = PhotoPairFactory::Get().GetPhotoAngleEnumFromString(photoangle_str);

    std::string annihilation_str = config.value("annihilation", "none");
    annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(annihilation_str);
    annihilation_def.multiplier = config.value("annihilation_multiplier", 1.0);
}


bool Utility::Definition::operator==(
    const Utility::Definition& utility_def) const {
    if (brems_def != utility_def.brems_def)
        return false;
    else if (compton_def != utility_def.compton_def)
        return false;
    else if (photo_def != utility_def.photo_def)
        return false;
    else if (epair_def != utility_def.epair_def)
        return false;
    else if (ioniz_def != utility_def.ioniz_def)
        return false;
    else if (mupair_def != utility_def.mupair_def)
        return false;
    else if (weak_def != utility_def.weak_def)
        return false;
    else if (photopair_def != utility_def.photopair_def)
        return false;
    else if (annihilation_def != utility_def.annihilation_def)
        return false;

    return true;
}

bool Utility::Definition::operator!=(
    const Utility::Definition& utility_def) const {
    return !(*this == utility_def);
}

Utility::Definition::~Definition() {}


// ------------------------------------------------------------------------- //
// OStream
// ------------------------------------------------------------------------- //

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, Utility const& utility)
{
    std::stringstream ss;
    ss << " Propagation Utility (" << &utility << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Particle Def:\n" << utility.particle_def_ << std::endl;
    os << "Medium:\n" << *utility.medium_ << std::endl;
    os << "EnergyCutSettings:\n" << utility.cut_settings_ << std::endl;

    os << Helper::Centered(60, "");
    return os;
}

} // namespace PROPOSAL

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

Utility::Utility(const ParticleDef& particle_def,
                 std::shared_ptr<const Medium> medium,
                 const EnergyCutSettings& cut_settings,
                 Definition utility_def)
    : particle_def_(particle_def)
    , medium_(medium)
    , cut_settings_(cut_settings)
    , crosssections_()
{
    if(utility_def.brems_def.parametrization!=BremsstrahlungFactory::Enum::None) {
        crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
                particle_def_, medium_, cut_settings_, utility_def.brems_def));
    }

    if(utility_def.photo_def.parametrization!=PhotonuclearFactory::Enum::None) {
        crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(
                particle_def_, medium_, cut_settings_,utility_def.photo_def));
    }

    if(utility_def.epair_def.parametrization!=EpairProductionFactory::Enum::None) {
        crosssections_.push_back(EpairProductionFactory::Get().CreateEpairProduction(
                particle_def_, medium_, cut_settings_, utility_def.epair_def));
    }

    if(utility_def.ioniz_def.parametrization!=IonizationFactory::Enum::None) {
        crosssections_.push_back(IonizationFactory::Get().CreateIonization(
                particle_def_, medium_, cut_settings_, utility_def.ioniz_def));
    }
    else{
        log_debug("No Ionization cross section chosen. For lepton propagation,Initialization may fail because no cross"
                  "section for small energies are available. You may have to enable Ionization or set a higher e_low"
                  "parameter for the particle.");
    }

    if(utility_def.annihilation_def.parametrization!=AnnihilationFactory::Enum::None) {
        crosssections_.push_back(AnnihilationFactory::Get().CreateAnnihilation(
                particle_def_, medium_, utility_def.annihilation_def));
        log_debug("Annihilation enabled");
    }

    if(utility_def.mupair_def.parametrization!=MupairProductionFactory::Enum::None) {
        crosssections_.push_back(MupairProductionFactory::Get().CreateMupairProduction(
            particle_def_, medium_, cut_settings_, utility_def.mupair_def));
        log_debug("Mupair Production enabled");
    }

    if(utility_def.weak_def.parametrization!=WeakInteractionFactory::Enum::None) {
        crosssections_.push_back(WeakInteractionFactory::Get().CreateWeakInteraction(
                    particle_def_, medium_, utility_def.weak_def));
        log_debug("Weak Interaction enabled");
    }

    // Photon interactions

    if(utility_def.compton_def.parametrization!=ComptonFactory::Enum::None) {
        crosssections_.push_back(ComptonFactory::Get().CreateCompton(
                particle_def_, medium_, cut_settings_, utility_def.compton_def));
        log_debug("Compton enabled");
    }

    if(utility_def.photopair_def.parametrization!=PhotoPairFactory::Enum::None) {
        crosssections_.push_back(PhotoPairFactory::Get().CreatePhotoPair(
                particle_def_, medium_, utility_def.photopair_def));
        log_debug("PhotoPairProduction enabled");
    }
}

Utility::Utility(const ParticleDef& particle_def,
                 std::shared_ptr<const Medium> medium,
                 const EnergyCutSettings& cut_settings,
                 Definition utility_def,
                 const InterpolationDef& interpolation_def)
    : particle_def_(particle_def)
    , medium_(medium)
    , cut_settings_(cut_settings)
    , crosssections_()
{
    if(utility_def.brems_def.parametrization!=BremsstrahlungFactory::Enum::None) {
        crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
                particle_def_, medium_, cut_settings_, utility_def.brems_def, interpolation_def));
    }

    if(utility_def.photo_def.parametrization!=PhotonuclearFactory::Enum::None) {
        crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(
                particle_def_, medium_, cut_settings_, utility_def.photo_def, interpolation_def));
    }

    if(utility_def.epair_def.parametrization!=EpairProductionFactory::Enum::None) {
        crosssections_.push_back(EpairProductionFactory::Get().CreateEpairProduction(
                particle_def_, medium_, cut_settings_, utility_def.epair_def, interpolation_def));
    }

    if(utility_def.ioniz_def.parametrization!=IonizationFactory::Enum::None) {
        crosssections_.push_back(IonizationFactory::Get().CreateIonization(
                particle_def_, medium_, cut_settings_, utility_def.ioniz_def, interpolation_def));
    }else{
        log_debug("No Ionization cross section chosen. For lepton propagation,Initialization may fail because no cross"
                  "section for small energies are available. You may have to enable Ionization or set a higher e_low"
                  "parameter for the particle.");
    }

    if(utility_def.annihilation_def.parametrization!=AnnihilationFactory::Enum::None) {
        crosssections_.push_back(AnnihilationFactory::Get().CreateAnnihilation(
                particle_def_, medium_, utility_def.annihilation_def, interpolation_def));
        log_debug("Annihilation enabled");
    }

    if(utility_def.mupair_def.parametrization!=MupairProductionFactory::Enum::None) {
        crosssections_.push_back(MupairProductionFactory::Get().CreateMupairProduction(
                    particle_def_, medium_, cut_settings_, utility_def.mupair_def, interpolation_def));
        log_debug("Mupair Production enabled");
    }

    if(utility_def.weak_def.parametrization!=WeakInteractionFactory::Enum::None) {
        crosssections_.push_back(WeakInteractionFactory::Get().CreateWeakInteraction(
                    particle_def_, medium_, utility_def.weak_def, interpolation_def));
        log_debug("Weak Interaction enabled");
    }

    // Photon interactions

    if(utility_def.compton_def.parametrization!=ComptonFactory::Enum::None) {
        crosssections_.push_back(ComptonFactory::Get().CreateCompton(
                particle_def_, medium_, cut_settings_, utility_def.compton_def, interpolation_def));
        log_debug("Compton enabled");
    }

    if(utility_def.photopair_def.parametrization!=PhotoPairFactory::Enum::None) {
        crosssections_.push_back(PhotoPairFactory::Get().CreatePhotoPair(
                particle_def_, medium_, utility_def.photopair_def, interpolation_def));
        log_debug("PhotoPairProduction enabled");
    }
}

Utility::Utility(const std::vector<CrossSection*>& crosssections) try
    : particle_def_(crosssections.at(0)->GetParametrization().GetParticleDef()),
      medium_(crosssections.at(0)->GetParametrization().GetMedium()),
      cut_settings_(crosssections.at(0)->GetParametrization().GetEnergyCuts()) {
    for (std::vector<CrossSection*>::const_iterator it = crosssections.begin();
         it != crosssections.end(); ++it) {
        if ((*it)->GetParametrization().GetParticleDef() != particle_def_) {
            log_fatal(
                "Particle definition of the cross section must be equal!");
        }

        if (*(*it)->GetParametrization().GetMedium() != *medium_) {
            log_fatal("Medium of the cross section must be equal!");
        }

        if ((*it)->GetParametrization().GetEnergyCuts() != cut_settings_) {
            log_fatal("Energy cuts of the cross section must be equal!");
        }

        crosssections_.push_back((*it)->clone());
    }
} catch (const std::out_of_range& e) {
    log_fatal("At least one cross section is needed for initializing utility.");
}

Utility::Utility(const Utility& collection)
    : particle_def_(collection.particle_def_),
      medium_(collection.medium_),
      cut_settings_(collection.cut_settings_),
      crosssections_(collection.crosssections_.size(), NULL) {
    for (unsigned int i = 0; i < crosssections_.size(); ++i) {
        crosssections_[i] = collection.crosssections_[i]->clone();
    }
}

Utility::~Utility() {

    for (std::vector<CrossSection*>::const_iterator iter =
             crosssections_.begin();
         iter != crosssections_.end(); ++iter) {
        delete *iter;
    }

    crosssections_.clear();
}

bool Utility::operator==(const Utility& utility) const {
    if (particle_def_ != utility.particle_def_)
        return false;
    else if (*medium_ != *utility.medium_)
        return false;
    else if (cut_settings_ != utility.cut_settings_)
        return false;
    else if (crosssections_.size() != utility.crosssections_.size())
        return false;

    for (unsigned int i = 0; i < crosssections_.size(); ++i) {
        if (*crosssections_[i] != *utility.crosssections_[i])
            return false;
    }

    return true;
}

bool Utility::operator!=(const Utility& utility) const {
    return !(*this == utility);
}

CrossSection* Utility::GetCrosssection(int typeId) const {
    for (auto& i : crosssections_) {
        if (i->GetTypeId() == typeId) {
            return i;
        }
    }
    log_fatal("No CrossSection found");
    return nullptr;
}


std::pair<double, int> Utility::StochasticLoss(
    double particle_energy, double rnd1, double rnd2, double rnd3)
{
    double total_rate = 0;
    double total_rate_weighted = 0;
    double rates_sum = 0;
    std::vector<double> rates;
    rates.resize(crosssections_.size());

    // return 0 and unknown, if there is no interaction
    std::pair<double, int> energy_loss;
    energy_loss.first = 0.;
    energy_loss.second = 0;

    for (unsigned int i = 0; i < crosssections_.size(); i++) {
        rates[i] = crosssections_[i]->CalculatedNdx(particle_energy, rnd2);
        total_rate += rates[i];
    }

    total_rate_weighted = total_rate * rnd1;

    log_debug("Total rate = %f, total rate weighted = %f", total_rate,
              total_rate_weighted);

    for (unsigned int i = 0; i < rates.size(); i++) {
        rates_sum += rates[i];

        if (rates_sum >= total_rate_weighted) {
            energy_loss.first = crosssections_[i]->CalculateStochasticLoss(
                particle_energy, rnd2, rnd3);
            energy_loss.second = crosssections_[i]->GetTypeId();

            break;
        }
    }

    return energy_loss;
}


/******************************************************************************
 *                            Utility Decorator                            *
 ******************************************************************************/

UtilityDecorator::UtilityDecorator(const Utility& utility)
    : utility_(*utility.clone()) {}

UtilityDecorator::UtilityDecorator(const UtilityDecorator& decorator)
    : utility_(*decorator.utility_.clone()) {}

UtilityDecorator::~UtilityDecorator() {}

bool UtilityDecorator::operator==(
    const UtilityDecorator& utility_decorator) const {
    if (typeid(*this) != typeid(utility_decorator))
        return false;
    if (utility_ != utility_decorator.utility_)
        return false;
    else
        return this->compare(utility_decorator);
}

bool UtilityDecorator::operator!=(
    const UtilityDecorator& utility_decorator) const {
    return !(*this == utility_decorator);
}

// ------------------------------------------------------------------------- //
double UtilityDecorator::FunctionToIntegral(double energy) {
    double result = 0.0;

    const std::vector<CrossSection*> crosssections =
        utility_.GetCrosssections();

    for (std::vector<CrossSection*>::const_iterator iter =
             crosssections.begin();
         iter != crosssections.end(); ++iter) {
        result += (*iter)->CalculatedEdx(energy);
    }

    return -1.0 / result;
}

double UtilityDecorator::Calculate(double ei,
                                   double ef,
                                   double rnd,
                                   const Vector3D& xi,
                                   const Vector3D& direction) {
    (void)xi;
    (void)direction;

    return this->Calculate(ei, ef, rnd);
}
