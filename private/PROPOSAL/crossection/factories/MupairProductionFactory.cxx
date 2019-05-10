
#include <algorithm>

#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/factories/MupairProductionFactory.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

MupairProductionFactory::MupairProductionFactory()
    : mupair_map_str_()
    , mupair_map_enum_()
    , string_enum_()
{
    Register("mupairkelnerkokoulinpetrukhin", KelnerKokoulinPetrukhin, std::make_pair(&MupairKelnerKokoulinPetrukhin::create, &MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>::create));
}

MupairProductionFactory::~MupairProductionFactory()
{
    string_enum_.clear();
    mupair_map_str_.clear();
    mupair_map_enum_.clear();
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* MupairProductionFactory::CreateMupairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def) const
{
    MupairMapEnum::const_iterator it = mupair_map_enum_.find(def.parametrization);

    if (it != mupair_map_enum_.end())
    {
        return new MupairIntegral(*it->second.first(particle_def, medium, cuts, def.multiplier));
    } else
    {
        log_fatal("MupairProduction %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
CrossSection* MupairProductionFactory::CreateMupairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    MupairMapEnum::const_iterator it = mupair_map_enum_.find(def.parametrization);

    if (it != mupair_map_enum_.end())
    {
        return new MupairInterpolant(*it->second.second(particle_def, medium, cuts, def.multiplier, interpolation_def), interpolation_def);
    } else
    {
        log_fatal("MupairProduction %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
void MupairProductionFactory::Register(const std::string& name,
                                     Enum enum_t,
                                     std::pair<RegisterFunction, RegisterFunctionInterpolant> create)
{
    mupair_map_str_[name]    = create;
    mupair_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// ------------------------------------------------------------------------- //
MupairProductionFactory::Enum MupairProductionFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("MupairProduction %s not registerd!", name.c_str());
        return None; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
std::string MupairProductionFactory::GetStringFromEnum(const MupairProductionFactory::Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("MupairProduction %s not registerd!", typeid(enum_t).name());
        return ""; // Just to prevent warinngs
    }
}
