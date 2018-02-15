
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include <algorithm>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringNoScattering.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

ScatteringFactory::ScatteringFactory()
    : registerd_enum()
    , registerd_str()
    , string_enum_()
{
    Register("moliere", Moliere);
    Register("highland", Highland);
    Register("highlandintegral", HighlandIntegral);
    Register("noscattering", NoScattering);
}

ScatteringFactory::~ScatteringFactory()
{
    // registerd_enum.clear();
    // registerd_str.clear();
    string_enum_.clear();
}

void ScatteringFactory::Register(const std::string& name, const Enum model)
{
    registerd_str.push_back(name);
    registerd_enum.push_back(model);
    string_enum_.insert(BimapStringEnum::value_type(name, model));
    // map_string_to_enum[name] = model;
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const std::string& name, Particle& particle, const Utility& utility, const InterpolationDef& interpolation_def)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);
    std::cout << name_lower << std::endl;

    std::vector<std::string>::const_iterator iter;
    iter = std::find(registerd_str.begin(), registerd_str.end(), name_lower);

    if (iter != registerd_str.end())
    {
        if (*iter == "highlandintegral")
        {
            return new ScatteringHighlandIntegral(particle, utility, interpolation_def);
        }
        else if (*iter == "moliere")
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        }
        else if (*iter == "highland")
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        }
        else if (*iter == "noscattering")
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        }
        else
        {
            log_fatal("Scattering %s not registerd!", name.c_str());
        }
    }
    else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const Enum model, Particle& particle, const Utility& utility, const InterpolationDef& interpolation_def)
{
    std::vector<Enum>::const_iterator iter;
    iter = std::find(registerd_enum.begin(), registerd_enum.end(), model);

    if (iter != registerd_enum.end())
    {
        if (*iter == HighlandIntegral)
        {
            return new ScatteringHighlandIntegral(particle, utility, interpolation_def);
        }
        else if (*iter == Moliere)
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        }
        else if (*iter == Highland)
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        }
        else if (*iter == NoScattering)
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        }
        else
        {
            log_fatal("Scattering %s not registerd!", typeid(model).name());
        }
    }
    else
    {
        log_fatal("Scattering %s not registerd!", typeid(model).name());
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const std::string& name, Particle& particle, const Utility& utility)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    std::vector<std::string>::const_iterator iter;
    iter = std::find(registerd_str.begin(), registerd_str.end(), name_lower);

    if (iter != registerd_str.end())
    {
        if (*iter == "highlandintegral")
        {
            return new ScatteringHighlandIntegral(particle, utility);
        }
        else if (*iter == "moliere")
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        }
        else if (*iter == "highland")
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        }
        else if (*iter == "noscattering")
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        }
        else
        {
            log_fatal("Scattering %s not registerd!", name.c_str());
        }
    }
    else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const Enum model, Particle& particle, const Utility& utility)
{
    std::vector<Enum>::const_iterator iter;
    iter = std::find(registerd_enum.begin(), registerd_enum.end(), model);

    if (iter != registerd_enum.end())
    {
        if (*iter == HighlandIntegral)
        {
            return new ScatteringHighlandIntegral(particle, utility);
        }
        else if (*iter == Moliere)
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        }
        else if (*iter == Highland)
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        }
        else if (*iter == NoScattering)
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        }
        else
        {
            log_fatal("Scattering %s not registerd!", typeid(model).name());
        }
    }
    else
    {
        log_fatal("Scattering %s not registerd!", typeid(model).name());
    }
}

// ------------------------------------------------------------------------- //
ScatteringFactory::Enum ScatteringFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Scattering model %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
std::string ScatteringFactory::GetStringFromEnum(const Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Scattering model %s not registerd!", typeid(enum_t).name());
    }
}

// ScatteringFactory::ScatteringFactory()
// {
//     // Register all media in lower case!
//
//     RegisterUtility("default", Default, &ScatteringDefault::create, &ScatteringDefault::create);
//     Register("moliere", Moliere, &ScatteringMoliere::create);
//     Register("moliere_first_order", MoliereFirstOrder, &ScatteringFirstOrder::create);
// }
//
// ScatteringFactory::~ScatteringFactory()
// {
//     scattering_map_str_.clear();
//     scattering_map_enum_.clear();
//
//     scattering_map_utility_str_.clear();
//     scattering_map_utility_enum_.clear();
//
//     map_string_to_enum.clear();
// }
//
// void ScatteringFactory::Register(const std::string& name, Enum model, RegisterFunction create)
// {
//     scattering_map_str_[name] = create;
//     scattering_map_enum_[model] = create;
//
//     map_string_to_enum[name] = model;
// }
//
// void ScatteringFactory::RegisterUtility(const std::string& name, Enum model, std::pair<RegisterFunctionUtility, RegisterFunctionUtilityInterpolant> create)
// {
//     scattering_map_utility_str_[name] = create;
//     scattering_map_utility_enum_[model] = create;
//
//     map_string_to_enum[name] = model;
// }

// ------------------------------------------------------------------------- //
// Creator functions
// ------------------------------------------------------------------------- //

// // ------------------------------------------------------------------------- //
// Scattering* ScatteringFactory::CreateScattering(const std::string& name, Particle& particle, Utility& utility)
// {
//     std::string name_lower = boost::algorithm::to_lower_copy(name);
//
//     ScatteringMapString::iterator it = scattering_map_str_.find(name_lower);
//     ScatteringMapUtiltiyString::iterator it_utility = scattering_map_utility_str_.find(name_lower);
//
//     if (it != scattering_map_str_.end())
//     {
//         return it->second(particle, utility.GetMedium());
//     }
//     else if (it_utility != scattering_map_utility_str_.end())
//     {
//         return it_utility->second(particle, utility);
//     }
//     else
//     {
//         log_fatal("Scattering %s not registerd!", name.c_str());
//     }
// }
//
// Scattering* ScatteringFactory::CreateScattering(const std::string& name, Particle& particle, Utility& utility, const InterpolationDef& interpolation_def)
// {
//     std::string name_lower = boost::algorithm::to_lower_copy(name);
//
//     ScatteringMapString::iterator it = scattering_map_str_.find(name_lower);
//     ScatteringMapUtiltiyString::iterator it_utility = scattering_map_utility_str_.find(name_lower);
//
//     if (it != scattering_map_str_.end())
//     {
//         return it->second(particle, utility.GetMedium());
//     }
//     else if (it_utility != scattering_map_utility_str_.end())
//     {
//         return it_utility->second(particle, utility, interpolation_def);
//     }
//     else
//     {
//         log_fatal("Scattering %s not registerd!", name.c_str());
//     }
// }
//
// Scattering* ScatteringFactory::CreateScattering(Enum model, Particle& particle, Utility& utility)
// {
//     ScatteringMapEnum::iterator it = scattering_map_enum_.find(model);
//     ScatteringMapUtiltiyEnum::iterator it_utility = scattering_map_utility_enum_.find(model);
//
//     if (it != scattering_map_enum_.end())
//     {
//         return it->second(particle, utility.GetMedium());
//     }
//     else if (it_utility != scattering_map_utility_enum_.end())
//     {
//         return it_utility->second(particle, utility);
//     }
//     else
//     {
//         log_fatal("Scattering %s not registerd!", typeid(model).name());
//     }
// }
