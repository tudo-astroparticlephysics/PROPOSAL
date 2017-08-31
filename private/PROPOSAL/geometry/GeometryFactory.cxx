
#include <boost/algorithm/string.hpp>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"

using namespace PROPOSAL;

GeometryFactory::GeometryFactory()
{
    // Geometry* (*create_sphere) () = &Sphere::create;
    // Geometry* (*create_box) () = &Box::create;
    // Geometry* (*create_cylinder) () = &Cylinder::create;
    //
    // Geometry* (*create_sphere_ptree) (boost::property_tree::ptree const&) = &Sphere::create;
    // Geometry* (*create_box_ptree) (boost::property_tree::ptree const&) = &Box::create;
    // Geometry* (*create_cylinder_ptree) (boost::property_tree::ptree const&) = &Cylinder::create;

    Register("sphere", &Sphere::create);
    Register("box", &Box::create);
    Register("cylinder", &Cylinder::create);

    // Register(GeometryModel::Default, &GeometryDefault::create);
    // Register(GeometryModel::Moliere, &GeometryMoliere::create);
    // Register(GeometryModel::MoliereFirstOrder, &GeometryFirstOrder::create);
}

GeometryFactory::~GeometryFactory()
{
    geometry_map.clear();
}

void GeometryFactory::Register(const std::string& name, RegisterFunction create)
{
    geometry_map[name] = create;
}

// void GeometryFactory::Register(GeometryModel::Enum model, RegisterFunction create)
// {
//     geometry_map[model] = create;
// }

Geometry* GeometryFactory::CreateGeometry(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    GeometryMap::iterator it = geometry_map.find(name_lower);

    if (it != geometry_map.end())
    {
        return it->second();
    } else
    {
        log_fatal("Geometry %s not registerd!", name.c_str());
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

    double x, y, z;
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
                break;
        }
        ++i;
    }

    Vector3D vec(x, y, z);

    // --------------------------------------------------------------------- //
    // Check type of geometry and set members
    // --------------------------------------------------------------------- //

    if (Sphere* sphere = dynamic_cast<Sphere*>(geometry))
    {
        double radius = pt.get<double>("radius");
        double inner_radius = pt.get<double>("inner_radius");

        sphere->SetPosition(vec);
        sphere->SetRadius(radius);
        sphere->SetInnerRadius(inner_radius);

        return sphere;
    }
    else if (Box* box = dynamic_cast<Box*>(geometry))
    {
        double x = pt.get<double>("x");
        double y = pt.get<double>("y");
        double z = pt.get<double>("z");

        box->SetPosition(vec);
        box->SetX(x);
        box->SetY(y);
        box->SetZ(z);

        return box;
    }
    else if (Cylinder* cylinder = dynamic_cast<Cylinder*>(geometry))
    {
        double radius = pt.get<double>("radius");
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
