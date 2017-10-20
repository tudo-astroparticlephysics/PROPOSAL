
#include <boost/algorithm/string.hpp>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"

using namespace PROPOSAL;

GeometryFactory::GeometryFactory()
{
    Register("sphere", Sphere, &Sphere::create);
    Register("box", Box, &Box::create);
    Register("cylinder", Cylinder, &Cylinder::create);
}

GeometryFactory::~GeometryFactory()
{
    geometry_map_str.clear();
    geometry_map_enum.clear();
}

void GeometryFactory::Register(const std::string& name, const Enum& num, RegisterFunction create)
{
    geometry_map_str[name] = create;
    geometry_map_enum[num] = create;
}

// void GeometryFactory::Register(GeometryModel::Enum model, RegisterFunction create)
// {
//     geometry_map[model] = create;
// }

Geometry* GeometryFactory::CreateGeometry(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    GeometryMapString::iterator it = geometry_map_str.find(name_lower);

    if (it != geometry_map_str.end())
    {
        return it->second();
    } else
    {
        log_fatal("Geometry %s not registerd!", name.c_str());
    }
}

Geometry* GeometryFactory::CreateGeometry(const Enum& num)
{
    GeometryMapEnum::iterator it = geometry_map_enum.find(num);

    if (it != geometry_map_enum.end())
    {
        return it->second();
    } else
    {
        log_fatal("Medium %s not registerd!", typeid(num).name());
    }
}

Geometry* GeometryFactory::CreateGeometry(const Definition& def)
{
    Geometry* geometry = CreateGeometry(def.shape);

    if (PROPOSAL::Sphere* sphere = dynamic_cast<PROPOSAL::Sphere*>(geometry))
    {
        sphere->SetPosition(def.position);
        sphere->SetRadius(def.radius);
        sphere->SetInnerRadius(def.inner_radius);
    }
    else if (PROPOSAL::Box* box = dynamic_cast<PROPOSAL::Box*>(geometry))
    {
        box->SetPosition(def.position);
        box->SetX(def.width);
        box->SetY(def.height);
        box->SetZ(def.depth);
    }

    else if (PROPOSAL::Cylinder* cylinder = dynamic_cast<PROPOSAL::Cylinder*>(geometry))
    {
        cylinder->SetPosition(def.position);
        cylinder->SetRadius(def.radius);
        cylinder->SetInnerRadius(def.inner_radius);
    }
    else
    {
        log_fatal("Geometry %s not registerd!", typeid(geometry).name());
    }

    return geometry;
}
