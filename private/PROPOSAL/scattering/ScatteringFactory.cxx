
#include <algorithm>

#include "PROPOSAL/Logging.h"
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
    string_enum_.insert(name, model);
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const std::string& name,
                                                const ParticleDef& particle_def,
                                                const Utility& utility,
                                                const InterpolationDef& interpolation_def)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    std::vector<std::string>::const_iterator iter;
    iter = std::find(registerd_str.begin(), registerd_str.end(), name_lower);

    if (iter != registerd_str.end())
    {
        if (*iter == "highlandintegral")
        {
            return new ScatteringHighlandIntegral(particle_def, utility, interpolation_def);
        } else if (*iter == "moliere")
        {
            return new ScatteringMoliere(particle_def, utility.GetMedium());
        } else if (*iter == "highland")
        {
            return new ScatteringHighland(particle_def, utility.GetMedium());
        } else if (*iter == "noscattering")
        {
            return new ScatteringNoScattering(particle_def, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", name.c_str());
            return NULL; // just to prevent warnings
        }
    } else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
        return NULL; // just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const Enum model,
                                                const ParticleDef& particle_def,
                                                const Utility& utility,
                                                const InterpolationDef& interpolation_def)
{
    std::vector<Enum>::const_iterator iter;
    iter = std::find(registerd_enum.begin(), registerd_enum.end(), model);

    if (iter != registerd_enum.end())
    {
        if (*iter == HighlandIntegral)
        {
            return new ScatteringHighlandIntegral(particle_def, utility, interpolation_def);
        } else if (*iter == Moliere)
        {
            return new ScatteringMoliere(particle_def, utility.GetMedium());
        } else if (*iter == Highland)
        {
            return new ScatteringHighland(particle_def, utility.GetMedium());
        } else if (*iter == NoScattering)
        {
            return new ScatteringNoScattering(particle_def, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", typeid(model).name());
            return NULL; // just to prevent warnings
        }
    } else
    {
        log_fatal("Scattering %s not registerd!", typeid(model).name());
        return NULL; // just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const std::string& name, const ParticleDef& particle_def, const Utility& utility)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    std::vector<std::string>::const_iterator iter;
    iter = std::find(registerd_str.begin(), registerd_str.end(), name_lower);

    if (iter != registerd_str.end())
    {
        if (*iter == "highlandintegral")
        {
            return new ScatteringHighlandIntegral(particle_def, utility);
        } else if (*iter == "moliere")
        {
            return new ScatteringMoliere(particle_def, utility.GetMedium());
        } else if (*iter == "highland")
        {
            return new ScatteringHighland(particle_def, utility.GetMedium());
        } else if (*iter == "noscattering")
        {
            return new ScatteringNoScattering(particle_def, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", name.c_str());
            return NULL; // just to prevent warnings
        }
    } else
    {
        log_fatal("Scattering %s not registerd!", name.c_str());
        return NULL; // just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
Scattering* ScatteringFactory::CreateScattering(const Enum model, const ParticleDef& particle_def, const Utility& utility)
{
    std::vector<Enum>::const_iterator iter;
    iter = std::find(registerd_enum.begin(), registerd_enum.end(), model);

    if (iter != registerd_enum.end())
    {
        if (*iter == HighlandIntegral)
        {
            return new ScatteringHighlandIntegral(particle_def, utility);
        } else if (*iter == Moliere)
        {
            return new ScatteringMoliere(particle_def, utility.GetMedium());
        } else if (*iter == Highland)
        {
            return new ScatteringHighland(particle_def, utility.GetMedium());
        } else if (*iter == NoScattering)
        {
            return new ScatteringNoScattering(particle_def, utility.GetMedium());
        } else
        {
            log_fatal("Scattering %s not registerd!", typeid(model).name());
            return NULL; // just to prevent warnings
        }
    } else
    {
        log_fatal("Scattering %s not registerd!", typeid(model).name());
        return NULL; // just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
ScatteringFactory::Enum ScatteringFactory::GetEnumFromString(const std::string& name)
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
        log_fatal("Scattering model %s not registerd!", name.c_str());
        return ScatteringFactory::None; // just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
std::string ScatteringFactory::GetStringFromEnum(const Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Scattering model %s not registerd!", typeid(enum_t).name());
        return ""; // just to prevent warnings
    }
}
