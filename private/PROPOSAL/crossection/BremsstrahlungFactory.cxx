
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/BremsstrahlungFactory.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

BremsstrahlungFactory::BremsstrahlungFactory()
{
    // Register all bremsstrahlung parametrizations in lower case!

    Register("PetrukhinShestakov", PetrukhinShestakov, &BremsPetrukhinShestakov::create);
    Register("KelnerKokoulinPetrukhin", KelnerKokoulinPetrukhin, &BremsKelnerKokoulinPetrukhin::create);
    Register("CompleteScreening", CompleteScreening, &BremsCompleteScreening::create);
    Register("AndreevBezrukovBugaev", AndreevBezrukovBugaev, &BremsAndreevBezrukovBugaev::create);
}

BremsstrahlungFactory::~BremsstrahlungFactory()
{
    bremsstrahlung_map_str_.clear();
    bremsstrahlung_map_enum_.clear();
    map_string_to_enum.clear();
}

void BremsstrahlungFactory::Register(const std::string& name, Enum enum_t, RegisterFunction create)
{
    bremsstrahlung_map_str_[name] = create;
    bremsstrahlung_map_enum_[enum_t] = create;
    map_string_to_enum[name] = enum_t;
}

// ------------------------------------------------------------------------- //
Parametrization* BremsstrahlungFactory::CreateParametrization(const std::string& name,
                                                              const ParticleDef& particle_def,
                                                              const Medium& medium,
                                                              const EnergyCutSettings& cuts,
                                                              Parametrization::Definition def) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BremsstrahlungMapString::const_iterator it = bremsstrahlung_map_str_.find(name_lower);

    if (it != bremsstrahlung_map_str_.end())
    {
        return it->second(particle_def, medium, cuts, def);
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Parametrization* BremsstrahlungFactory::CreateParametrization(Enum enum_t,
                                                              const ParticleDef& particle_def,
                                                              const Medium& medium,
                                                              const EnergyCutSettings& cuts,
                                                              Parametrization::Definition def) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(enum_t);

    if (it != bremsstrahlung_map_enum_.end())
    {
        return it->second(particle_def, medium, cuts, def);
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const std::string& name,
                                                              const ParticleDef& particle_def,
                                                              const Medium& medium,
                                                              const EnergyCutSettings& cuts,
                                                              Parametrization::Definition def,
                                                              bool interpolate) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BremsstrahlungMapString::const_iterator it = bremsstrahlung_map_str_.find(name_lower);

    if (it != bremsstrahlung_map_str_.end())
    {
        if (interpolate)
        {
            return new BremsInterpolant(*it->second(particle_def, medium, cuts, def));
        }
        else
        {
            return new BremsIntegral(*it->second(particle_def, medium, cuts, def));
        }
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(Enum enum_t,
                                                          const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          Parametrization::Definition def,
                                                          bool interpolate) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(enum_t);

    if (it != bremsstrahlung_map_enum_.end())
    {
        if (interpolate)
        {
            return new BremsInterpolant(*it->second(particle_def, medium, cuts, def));
        }
        else
        {
            return new BremsIntegral(*it->second(particle_def, medium, cuts, def));
        }
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
BremsstrahlungFactory::Enum BremsstrahlungFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    MapStringToEnum::iterator it = map_string_to_enum.find(name_lower);

    if (it != map_string_to_enum.end())
    {
        return it->second;
    } else
    {
        log_fatal("Bremsstrahlung %s not registerd!", name.c_str());
    }
}
