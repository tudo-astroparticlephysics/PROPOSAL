
#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>

#include "PROPOSAL/geometry/Geometry.h"

namespace PROPOSAL
{

class GeometryFactory
{
    public:

    typedef boost::function<Geometry* (void)> RegisterFunction;
    // typedef boost::function<Geometry* (boost::property_tree::ptree const&)> RegisterFunctionPTree;
    typedef std::map<std::string, RegisterFunction> GeometryMap;
    // typedef std::map<std::string, RegisterFunctionPTree> GeometryMapPTree;
    // typedef std::map<ScatteringModel::Enum, boost::function<Scattering* (void)> > ScatteringMapEnum;

    // void Register(const std::string& name, RegisterFunction);
    void Register(const std::string& name, RegisterFunction);

    Geometry* CreateGeometry(const std::string&);
    Geometry* CreateGeometry(boost::property_tree::ptree const& pt);

    static GeometryFactory& Get()
    {
        static GeometryFactory instance;
        return instance;
    }

    private:
    GeometryFactory();
    ~GeometryFactory();

    GeometryMap geometry_map;
    // GeometryMapPTree geometry_map_ptree;
    // ScatteringMapEnum scattering_map_enum_;
};

} /* PROPOSAL */

