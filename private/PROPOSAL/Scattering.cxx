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

void Scattering::Scatter(PROPOSALParticle& particle, const std::vector<CrossSections*>& cross_sections, double dr, double ei, double ef)
{
    double sz,tz;

    RandomAngles random_angles = CalculateRandomAngle(particle, cross_sections, dr, ei, ef);

    sz = sqrt(std::max(1.-(random_angles.sx*random_angles.sx+random_angles.sy*random_angles.sy), 0.));
    tz = sqrt(std::max(1.-(random_angles.tx*random_angles.tx+random_angles.ty*random_angles.ty), 0.));

    Vector3D position;
    Vector3D direction;

    long double sinth, costh,sinph,cosph;
    sinth = (long double) sin(particle.GetDirection().GetTheta());
    costh = (long double) cos(particle.GetDirection().GetTheta());
    sinph = (long double) sin(particle.GetDirection().GetPhi());
    cosph = (long double) cos(particle.GetDirection().GetPhi());

    position = particle.GetPosition();

    // Rotation towards all tree axes
    direction = sz*particle.GetDirection();
    direction = direction + random_angles.sx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + random_angles.sy*Vector3D(-sinph, cosph, 0.);

    position = position + dr*direction;

    // Rotation towards all tree axes
    direction = tz*particle.GetDirection();
    direction = direction + random_angles.tx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + random_angles.ty*Vector3D(-sinph, cosph, 0.);

    direction.CalculateSphericalCoordinates();

    particle.SetPosition(position);
    particle.SetDirection(direction);
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
