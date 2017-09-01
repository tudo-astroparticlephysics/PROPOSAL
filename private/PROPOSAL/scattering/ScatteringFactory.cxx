
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/Output.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/scattering/ScatteringFirstOrder.h"

using namespace PROPOSAL;

ScatteringFactory::ScatteringFactory()
{
    // Register all media in lower case!

    Register("default", &ScatteringDefault::create);
    Register("moliere", &ScatteringMoliere::create);
    Register("moliere_first_order", &ScatteringFirstOrder::create);

    Register("default", Default);
    Register("moliere", Moliere);
    Register("moliere_first_order", MoliereFirstOrder);

    Register(Default, &ScatteringDefault::create);
    Register(Moliere, &ScatteringMoliere::create);
    Register(MoliereFirstOrder, &ScatteringFirstOrder::create);
}

ScatteringFactory::~ScatteringFactory()
{
    scattering_map_str_.clear();
    scattering_map_enum_.clear();
    map_string_to_enum.clear();
}

void ScatteringFactory::Register(const std::string& name, RegisterFunction create)
{
    scattering_map_str_[name] = create;
}

void ScatteringFactory::Register(Enum model, RegisterFunction create)
{
    scattering_map_enum_[model] = create;
}

void ScatteringFactory::Register(const std::string& name, Enum model)
{
    map_string_to_enum[name] = model;
}

Scattering* ScatteringFactory::CreateScattering(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    ScatteringMapString::iterator it = scattering_map_str_.find(name_lower);

    if (it != scattering_map_str_.end())
    {
        return it->second();
    } else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
    }
}

Scattering* ScatteringFactory::CreateScattering(Enum model)
{
    ScatteringMapEnum::iterator it = scattering_map_enum_.find(model);

    if (it != scattering_map_enum_.end())
    {
        return it->second();
    } else
    {
        log_fatal("Scattering %s not registerd!", typeid(model).name());
    }
}

ScatteringFactory::Enum ScatteringFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);
    MapStringToEnum::iterator it = map_string_to_enum.find(name_lower);

    if (it != map_string_to_enum.end())
    {
        return it->second;
    } else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
    }
}
