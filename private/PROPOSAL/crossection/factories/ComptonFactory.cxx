
#include <algorithm>
#include <stdexcept>

#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/ComptonInterpolant.h"
#include "PROPOSAL/crossection/factories/ComptonFactory.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

ComptonFactory::ComptonFactory()
    : compton_map_str_()
    , compton_map_enum_()
    , string_enum_()
{
    // Register all compton parametrizations in lower case!

    Register("comptonkleinnishina", KleinNishina, &ComptonKleinNishina::create);
    Register("none", None, nullptr); // empty parametrization
}

// ------------------------------------------------------------------------- //
ComptonFactory::~ComptonFactory()
{
    compton_map_str_.clear();
    compton_map_enum_.clear();
    string_enum_.clear();
}

// ------------------------------------------------------------------------- //
void ComptonFactory::Register(
    const std::string& name, Enum enum_t, RegisterFunction create)
{
    compton_map_str_[name] = create;
    compton_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* ComptonFactory::CreateCompton(const ParticleDef& particle_def,
    std::shared_ptr<const Medium> medium,
    std::shared_ptr<const EnergyCutSettings> cuts, const Definition& def,
    std::shared_ptr<const InterpolationDef> interpolation_def = nullptr) const
{
    ComptonMapEnum::const_iterator it
        = compton_map_enum_.find(def.parametrization);

    if (it != compton_map_enum_.end()) {
        if (interpolation_def) {
            return new ComptonInterpolant(
                *it->second(particle_def, medium, def.multiplier), cuts,
                *interpolation_def);
        }
        return new ComptonIntegral(
            *it->second(particle_def, medium, def.multiplier), cuts);
    }
    std::invalid_argument("Compton not registered!");
}

CrossSection* ComptonFactory::CreateCompton(const Compton& parametrization,
    std::shared_ptr<const EnergyCutSettings> cuts,
    std::shared_ptr<const InterpolationDef> interpolation_def = nullptr) const
{
    if (interpolation_def) {
        return new ComptonInterpolant(parametrization, cuts, *interpolation_def);
    }
    return new ComptonIntegral(parametrization, cuts);
}

// ------------------------------------------------------------------------- //
ComptonFactory::Enum ComptonFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end()) {
        return it->second;
    } else {
        log_fatal("Compton %s not registered!", name.c_str());
        return ComptonFactory::Fail; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
std::string ComptonFactory::GetStringFromEnum(
    const ComptonFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end()) {
        return it->second;
    } else {
        log_fatal("Compton %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
