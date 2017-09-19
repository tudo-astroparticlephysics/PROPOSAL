
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"

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
CrossSection* BremsstrahlungFactory::CreateBremsstrahlungIntegral(const std::string& name,
                                                          const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          double multiplier,
                                                          bool lpm) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BremsstrahlungMapString::const_iterator it = bremsstrahlung_map_str_.find(name);

    if (it != bremsstrahlung_map_str_.end())
    {
        return new BremsIntegral(*it->second(particle_def, medium, cuts, multiplier, lpm));
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}


// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlungIntegral(Enum enum_t,
                                                          const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          double multiplier,
                                                          bool lpm) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(enum_t);

    if (it != bremsstrahlung_map_enum_.end())
    {
        return new BremsIntegral(*it->second(particle_def, medium, cuts, multiplier, lpm));
    }
    else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlungInterpolant(const std::string& name,
                                                          const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          double multiplier,
                                                          bool lpm,
                                                          InterpolationDef def) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BremsstrahlungMapString::const_iterator it = bremsstrahlung_map_str_.find(name);

    if (it != bremsstrahlung_map_str_.end())
    {
        return new BremsInterpolant(*it->second(particle_def, medium, cuts, multiplier, lpm), def);
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}


// ------------------------------------------------------------------------- //
CrossSection* BremsstrahlungFactory::CreateBremsstrahlungInterpolant(Enum enum_t,
                                                          const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          double multiplier,
                                                          bool lpm,
                                                          InterpolationDef def) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(enum_t);

    if (it != bremsstrahlung_map_enum_.end())
    {
        return new BremsInterpolant(*it->second(particle_def, medium, cuts, multiplier, lpm), def);
    }
    else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(enum_t).name());
    }
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const Enum enum_t,
                                 const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 bool lpm,
                                 bool interpolate,
                                 InterpolationDef def) const
{
    BremsstrahlungMapEnum::const_iterator it = bremsstrahlung_map_enum_.find(enum_t);

    if (it != bremsstrahlung_map_enum_.end())
    {
        if (interpolate)
        {
            return new BremsInterpolant(*it->second(particle_def, medium, cuts, multiplier, lpm), def);
        }
        else
        {
            return new BremsIntegral(*it->second(particle_def, medium, cuts, multiplier, lpm));
        }
    }
    else
    {
        log_fatal("Bremsstrahlung %s not registerd!", typeid(enum_t).name());
    }
}

CrossSection* BremsstrahlungFactory::CreateBremsstrahlung(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          const Definition& def,
                                                          bool interpolate,
                                                          InterpolationDef interpolation_def) const
{
    return CreateBremsstrahlung(def.parametrization, particle_def, medium, cuts, def.multiplier, def.lpm_effect, interpolate, interpolation_def);
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
