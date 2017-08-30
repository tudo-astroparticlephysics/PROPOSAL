
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

MediumFactory::MediumFactory()
{
    // Register all media in lower case!

    Register("water", &Water::create);
    Register("ice", &Ice::create);
    Register("salt", &Salt::create);
    Register("standardrock", &StandardRock::create);
    Register("frejusrock", &FrejusRock::create);
    Register("iron", &Iron::create);
    Register("hydrogen", &Hydrogen::create);
    Register("lead", &Lead::create);
    Register("copper", &Copper::create);
    Register("uranium", &Uranium::create);
    Register("air", &Air::create);
    Register("paraffin", &AntaresWater::create);
}

void MediumFactory::Register(const std::string& name, Medium* (*create)(void))
{
    medium_map[name] = create;
}

Medium* MediumFactory::CreateMedium(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    std::map<std::string, Medium* (*)(void)>::iterator it = medium_map.find(name_lower);

    if (it != medium_map.end())
    {
        return it->second();
    } else
    {
        log_fatal("Medium %s not registerd!", name.c_str());
    }
}
