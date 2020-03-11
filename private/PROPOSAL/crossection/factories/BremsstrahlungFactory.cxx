
#include <algorithm>
#include <stdexcept>

#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

BremsstrahlungFactory::BremsstrahlungFactory()
    : bremsstrahlung_map_str_()
    , bremsstrahlung_map_enum_()
    , string_enum_()
{
    // Register all bremsstrahlung parametrizations in lower case!

    Register("bremspetrukhinshestakov", PetrukhinShestakov, &BremsPetrukhinShestakov::create);
    Register("bremskelnerkokoulinpetrukhin", KelnerKokoulinPetrukhin, &BremsKelnerKokoulinPetrukhin::create);
    Register("bremscompletescreening", CompleteScreening, &BremsCompleteScreening::create);
    Register("bremsandreevbezrukovbugaev", AndreevBezrukovBugaev, &BremsAndreevBezrukovBugaev::create);
    Register("bremssandrocksoedingreksorhode", SandrockSoedingreksoRhode, &BremsSandrockSoedingreksoRhode::create);
    Register("bremselectronscreening", ElectronScreening, &BremsElectronScreening::create);
    Register("none", None, nullptr); //empty parametrization
}

// ------------------------------------------------------------------------- //
BremsstrahlungFactory::~BremsstrahlungFactory()
{
    bremsstrahlung_map_str_.clear();
    bremsstrahlung_map_enum_.clear();
    string_enum_.clear();
}

// ------------------------------------------------------------------------- //
void BremsstrahlungFactory::Register(const std::string& name, Enum enum_t, RegisterFunction create)
{
    bremsstrahlung_map_str_[name]    = create;
    bremsstrahlung_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const ParticleDef& particle_def,
                                                          std::shared_ptr<const Medium> medium,
                                                          std::shared_ptr<const EnergyCutSettings> cuts,
                                                          const Definition& def,
                                                          std::shared_ptr<const InterpolationDef> interpolation_def = nullptr) const
{

    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(def.parametrization);

    if (it != bremsstrahlung_map_enum_.end())
    {
        if(interpolation_def)
            return new BremsInterpolant(*it->second(particle_def, medium, def.multiplier, def.lpm_effect), cuts,
                                        interpolation_def);

        return new BremsIntegral(*it->second(particle_def, medium, def.multiplier, def.lpm_effect), cuts);
    }

    std::invalid_argument("Bremsstrahlung not registered!");
}

CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const Bremsstrahlung& param,
                                                          std::shared_ptr<const EnergyCutSettings> cuts,
                                                          std::shared_ptr<const InterpolationDef> interpolation_def = nullptr) const
{
    if(interpolation_def){
        return new BremsInterpolant(param, cuts, interpolation_def);
    }

    return new BremsIntegral(param, cuts);
}

// ------------------------------------------------------------------------- //
BremsstrahlungFactory::Enum BremsstrahlungFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Bremsstrahlung %s not registered!", name.c_str());
        return BremsstrahlungFactory::Fail; // Just to prevent warinngs

    }
}

// ------------------------------------------------------------------------- //
std::string BremsstrahlungFactory::GetStringFromEnum(const BremsstrahlungFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Bremsstrahlung %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
