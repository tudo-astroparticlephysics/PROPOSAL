
#pragma once

#include <boost/function.hpp>
#include <map>

#include "PROPOSAL/geometry/Geometry.h"

namespace PROPOSAL {

class GeometryFactory
{
public:
    enum Enum
    {
        Sphere = 0,
        Box,
        Cylinder
    };

    struct Definition
    {
        Definition()
            : shape(Sphere)
            , position()
            , inner_radius(0.0)
            , radius(0.0)
            , width(0.0)
            , height(0.0)
            , depth(0.0)
        {
        }

        Enum shape;
        Vector3D position;
        double inner_radius;
        double radius;
        double width;
        double height;
        double depth;
    };

    typedef boost::function<Geometry*(void)> RegisterFunction;
    typedef std::map<std::string, RegisterFunction> GeometryMapString;
    typedef std::map<Enum, RegisterFunction> GeometryMapEnum;

    void Register(const std::string&, const Enum&, RegisterFunction);

    Geometry* CreateGeometry(const std::string&);
    Geometry* CreateGeometry(const Enum&);
    Geometry* CreateGeometry(const Definition&);

    static GeometryFactory& Get()
    {
        static GeometryFactory instance;
        return instance;
    }

private:
    GeometryFactory();
    ~GeometryFactory();

    GeometryMapString geometry_map_str;
    GeometryMapEnum geometry_map_enum;
};

} // namespace PROPOSAL
