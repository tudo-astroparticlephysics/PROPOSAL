#include <algorithm>
#include <stdexcept>

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

IonizationFactory::IonizationFactory()
    : ioniz_map_str_()
    , ioniz_map_enum_()
    , string_enum_()
{
    Register(
        "ionizbetheblochrossi", BetheBlochRossi, &IonizBetheBlochRossi::create);
    Register("ionizbergerseltzerbhabha", IonizBergerSeltzerBhabha,
        &IonizBergerSeltzerBhabha::create);
    Register("ionizbergerseltzermoller", IonizBergerSeltzerMoller,
        &IonizBergerSeltzerMoller::create);
}

// ------------------------------------------------------------------------- //
IonizationFactory::~IonizationFactory()
{
    ioniz_map_str_.clear();
    ioniz_map_enum_.clear();
    string_enum_.clear();
}

// ------------------------------------------------------------------------- //
void IonizationFactory::Register(
    const std::string& name, Enum enum_t, RegisterFunction create)
{
    ioniz_map_str_[name] = create;
    ioniz_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* IonizationFactory::CreateIonization(
    const ParticleDef& particle_def, std::shared_ptr<const Medium> medium,
    std::shared_ptr<const EnergyCutSettings> cuts, const Definition& def,
    std::shared_ptr<const InterpolationDef> interpolation_def) const
{
    IonizMapEnum::const_iterator it = ioniz_map_enum_.find(def.parametrization);

    if (it != ioniz_map_enum_.end()) {
        if (interpolation_def) {
            return new IonizInterpolant(
                *it->second(particle_def, medium, cuts, def.multiplier), cuts,
                *interpolation_def);
        }
        return new IonizIntegral(
            *it->second(particle_def, medium, cuts, def.multiplier), cuts);
    }
    std::invalid_argument("Ionization not registered!");
}

CrossSection* IonizationFactory::CreateIonization(
    const Ionization& parametrization,
    std::shared_ptr<const EnergyCutSettings> cuts,
    std::shared_ptr<const InterpolationDef> interpolation_def) const
{
    if (interpolation_def) {
        return new IonizInterpolant(parametrization, cuts, *interpolation_def);
    }
    return new IonizIntegral(parametrization, cuts);
}

// ------------------------------------------------------------------------- //
IonizationFactory::Enum IonizationFactory::GetEnumFromString(
    const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end()) {
        return it->second;
    } else {
        log_fatal("Ionization %s not registered!", name.c_str());
        return IonizationFactory::Fail; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
std::string IonizationFactory::GetStringFromEnum(
    const IonizationFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end()) {
        return it->second;
    } else {
        log_fatal("Ionization %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
