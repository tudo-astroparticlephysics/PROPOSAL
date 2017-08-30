/*! \file   Scattering.cxx
*   \brief  Source filefor the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
**/

#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/Output.h"
#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/ScatteringDefault.h"
#include "PROPOSAL/ScatteringMoliere.h"
#include "PROPOSAL/ScatteringFirstOrder.h"

using namespace PROPOSAL;

/******************************************************************************
*                                 Scattering                                  *
******************************************************************************/


Scattering::Scattering()
{
}

Scattering::~Scattering()
{
}

/******************************************************************************
*                             ScatteringFactory                               *
******************************************************************************/


ScatteringFactory::ScatteringFactory()
{
    // Register all media in lower case!

    Register("default", &ScatteringDefault::create);
    Register("moliere", &ScatteringMoliere::create);
    Register("moliere_first_order", &ScatteringFirstOrder::create);

    Register("default", ScatteringModel::Default);
    Register("moliere", ScatteringModel::Moliere);
    Register("moliere_first_order", ScatteringModel::MoliereFirstOrder);

    Register(ScatteringModel::Default, &ScatteringDefault::create);
    Register(ScatteringModel::Moliere, &ScatteringMoliere::create);
    Register(ScatteringModel::MoliereFirstOrder, &ScatteringFirstOrder::create);
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

void ScatteringFactory::Register(ScatteringModel::Enum model, RegisterFunction create)
{
    scattering_map_enum_[model] = create;
}

void ScatteringFactory::Register(const std::string& name, ScatteringModel::Enum model)
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

Scattering* ScatteringFactory::CreateScattering(ScatteringModel::Enum model)
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

ScatteringFactory::ScatteringModel::Enum ScatteringFactory::GetEnumFromString(const std::string& name)
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
