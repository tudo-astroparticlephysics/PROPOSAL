
#pragma once

#include <boost/function.hpp>

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/scattering/Scattering.h"

namespace  PROPOSAL
{

class ScatteringFactory
{
    public:

    struct ScatteringModel
    {
        enum Enum
        {
            Default = 0,
            Moliere,
            MoliereFirstOrder
        };
    };

    typedef boost::function<Scattering* (void)> RegisterFunction;
    typedef std::map<std::string, boost::function<Scattering* (void)> > ScatteringMapString;
    typedef std::map<ScatteringModel::Enum, boost::function<Scattering* (void)> > ScatteringMapEnum;
    typedef std::map<std::string, ScatteringModel::Enum> MapStringToEnum;

    Scattering* CreateScattering(const std::string&);
    Scattering* CreateScattering(const ScatteringModel::Enum);

    ScatteringModel::Enum GetEnumFromString(const std::string&);

    static ScatteringFactory& Get()
    {
        static ScatteringFactory instance;
        return instance;
    }

    private:
    ScatteringFactory();
    ~ScatteringFactory();

    void Register(const std::string& name, RegisterFunction);
    void Register(ScatteringModel::Enum, RegisterFunction);
    void Register(const std::string&, ScatteringModel::Enum);

    // std::map<std::string, Scattering* (*)(void)> scattering_map;
    ScatteringMapString scattering_map_str_;
    ScatteringMapEnum scattering_map_enum_;
    MapStringToEnum map_string_to_enum;
};

} /*  PROPOSAL */

