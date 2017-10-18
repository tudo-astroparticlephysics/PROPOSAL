
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

MediumFactory::MediumFactory()
    : medium_map_str()
    , medium_map_enum()
    , string_enum_()
{
    // Register all media in lower case!

    Register("water", Water, &Water::create);
    Register("ice", Ice, &Ice::create);
    Register("salt", Salt, &Salt::create);
    Register("standardrock", StandardRock, &StandardRock::create);
    Register("frejusrock", FrejusRock, &FrejusRock::create);
    Register("iron", Iron, &Iron::create);
    Register("hydrogen", Hydrogen, &Hydrogen::create);
    Register("lead", Lead, &Lead::create);
    Register("copper", Copper, &Copper::create);
    Register("uranium", Uranium, &Uranium::create);
    Register("air", Air, &Air::create);
    Register("paraffin", Paraffin, &AntaresWater::create);
}

MediumFactory::~MediumFactory()
{
    medium_map_str.clear();
    medium_map_enum.clear();
    string_enum_.clear();
}

void MediumFactory::Register(const std::string& name, const Enum& enum_t, RegisterFunction create)
{
    medium_map_str[name] = create;
    medium_map_enum[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// ------------------------------------------------------------------------- //
Medium* MediumFactory::CreateMedium(const std::string& name, double density_correction)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    MediumMapString::iterator it = medium_map_str.find(name_lower);

    if (it != medium_map_str.end())
    {
        return it->second(density_correction);
    } else
    {
        log_fatal("Medium %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Medium* MediumFactory::CreateMedium(const Enum& med, double density_correction)
{
    MediumMapEnum::iterator it = medium_map_enum.find(med);

    if (it != medium_map_enum.end())
    {
        return it->second(density_correction);
    } else
    {
        log_fatal("Medium %s not registerd!", typeid(med).name());
    }
}

// ------------------------------------------------------------------------- //
Medium* MediumFactory::CreateMedium(Definition def)
{
    MediumMapEnum::iterator it = medium_map_enum.find(def.type);

    if (it != medium_map_enum.end())
    {
        return it->second(def.density_correction);
    } else
    {
        log_fatal("Medium %s not registerd!", typeid(def.type).name());
    }
}

// ------------------------------------------------------------------------- //
MediumFactory::Enum MediumFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Medium %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
std::string MediumFactory::GetStringFromEnum(const MediumFactory::Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Medium %s not registerd!", typeid(enum_t).name());
    }
}
