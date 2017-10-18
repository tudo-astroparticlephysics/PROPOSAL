
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

Geometry* GeometryFactory::CreateGeometry(boost::property_tree::ptree const& pt)
{
    // --------------------------------------------------------------------- //
    // Get geometry from default constructor
    // --------------------------------------------------------------------- //

    std::string name = pt.get<std::string>("shape");
    Geometry* geometry = CreateGeometry(name);

    // --------------------------------------------------------------------- //
    // Get the position vector from the property tree
    // --------------------------------------------------------------------- //

    double x = 0;
    double y = 0;
    double z = 0;

    boost::property_tree::ptree child = pt.get_child("origin");

    if (child.size() != 3)
    {
        log_fatal("Cannot initialize Vector3D by Property Tree! Something wrong with the config file?");
    }

    int i = 0;
    for(boost::property_tree::ptree::const_iterator it = child.begin(); it != child.end(); ++it)
    {
        double coord = it->second.get_value<double>();

        switch (i)
        {
            case 0:
                x = coord;
                break;
            case 1:
                y = coord;
                break;
            case 2:
                z = coord;
                break;
            default:
                // Do nothing
                break;
        }
        ++i;
    }

    Vector3D vec(x, y, z);

    // --------------------------------------------------------------------- //
    // Check type of geometry and set members
    // --------------------------------------------------------------------- //

    if (PROPOSAL::Sphere* sphere = dynamic_cast<PROPOSAL::Sphere*>(geometry))
    {
        double radius = pt.get<double>("outer_radius");
        double inner_radius = pt.get<double>("inner_radius");

        sphere->SetPosition(vec);
        sphere->SetRadius(radius);
        sphere->SetInnerRadius(inner_radius);

        return sphere;
    }
    else if (PROPOSAL::Box* box = dynamic_cast<PROPOSAL::Box*>(geometry))
    {
        double x = pt.get<double>("lenght");
        double y = pt.get<double>("width");
        double z = pt.get<double>("height");

        box->SetPosition(vec);
        box->SetX(x);
        box->SetY(y);
        box->SetZ(z);

        return box;
    }
    else if (PROPOSAL::Cylinder* cylinder = dynamic_cast<PROPOSAL::Cylinder*>(geometry))
    {
        double radius = pt.get<double>("outer_radius");
        double inner_radius = pt.get<double>("inner_radius");
        double z = pt.get<double>("z");

        cylinder->SetPosition(vec);
        cylinder->SetRadius(radius);
        cylinder->SetInnerRadius(inner_radius);
        cylinder->SetZ(z);

        return cylinder;
    }
    else
    {
        log_fatal("GeometryFactory could not create %s with its properties!", name.c_str());
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
