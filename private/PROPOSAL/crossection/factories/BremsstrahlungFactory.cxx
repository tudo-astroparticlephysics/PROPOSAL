
#include <algorithm>

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
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          const Definition& def) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(def.parametrization);

    if (it != bremsstrahlung_map_enum_.end())
    {
        return new BremsIntegral(*it->second(particle_def, medium, cuts, def.multiplier, def.lpm_effect));
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          const Definition& def,
                                                          InterpolationDef interpolation_def) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(def.parametrization);

    if (it != bremsstrahlung_map_enum_.end())
    {
        return new BremsInterpolant(*it->second(particle_def, medium, cuts, def.multiplier, def.lpm_effect),
                                    interpolation_def);
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
BremsstrahlungFactory::Enum BremsstrahlungFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", name.c_str());
        return BremsstrahlungFactory::None; // Just to prevent warinngs

    }
}

// ------------------------------------------------------------------------- //
std::string BremsstrahlungFactory::GetStringFromEnum(const BremsstrahlungFactory::Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(enum_t).name());
        return ""; // Just to prevent warinngs
    }
}
