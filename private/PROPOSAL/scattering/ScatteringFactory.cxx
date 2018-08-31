
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include <algorithm>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
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
    string_enum_.clear();
}

void ScatteringFactory::Register(const std::string& name, const Enum model)
{
    registerd_str.push_back(name);
    registerd_enum.push_back(model);
    string_enum_.insert(BimapStringEnum::value_type(name, model));
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const std::string& name,
                                                Particle& particle,
                                                const Utility& utility,
                                                const InterpolationDef& interpolation_def)
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
        } else if (*iter == "moliere")
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        } else if (*iter == "highland")
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        } else if (*iter == "noscattering")
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", name.c_str());
        }
    } else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const Enum model,
                                                Particle& particle,
                                                const Utility& utility,
                                                const InterpolationDef& interpolation_def)
{
    std::vector<Enum>::const_iterator iter;
    iter = std::find(registerd_enum.begin(), registerd_enum.end(), model);

    if (iter != registerd_enum.end())
    {
        if (*iter == HighlandIntegral)
        {
            return new ScatteringHighlandIntegral(particle, utility, interpolation_def);
        } else if (*iter == Moliere)
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        } else if (*iter == Highland)
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        } else if (*iter == NoScattering)
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", typeid(model).name());
        }
    } else
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
        } else if (*iter == "moliere")
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        } else if (*iter == "highland")
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        } else if (*iter == "noscattering")
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", name.c_str());
        }
    } else
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
        } else if (*iter == Moliere)
        {
            return new ScatteringMoliere(particle, utility.GetMedium());
        } else if (*iter == Highland)
        {
            return new ScatteringHighland(particle, utility.GetMedium());
        } else if (*iter == NoScattering)
        {
            return new ScatteringNoScattering(particle, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", typeid(model).name());
        }
    } else
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
