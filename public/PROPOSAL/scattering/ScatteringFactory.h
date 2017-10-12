
#pragma once

// #include <boost/function.hpp>

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/scattering/Scattering.h"

namespace  PROPOSAL
{

class Medium;
class Utility;
struct InterpolationDef;

class ScatteringFactory
{
    public:

    enum Enum
    {
        Default = 0,
        Moliere,
        Highland,
        NoScattering
    };

    // typedef boost::function<Scattering* (Particle&, const Medium&)> RegisterFunction;
    // typedef boost::function<Scattering* (Particle&, Utility&)> RegisterFunctionUtility;
    // typedef boost::function<Scattering* (Particle&, Utility&)> RegisterFunctionUtilityInterpolant;
    //
    // typedef std::map<std::string, RegisterFunction> ScatteringMapString;
    // typedef std::map<Enum, RegisterFunction> ScatteringMapEnum;
    //
    // typedef std::map<std::string, std::pair<RegisterFunctionUtility, RegisterFunctionUtilityInterpolant> > ScatteringMapUtiltiyString;
    // typedef std::map<Enum, std::pair<RegisterFunctionUtility, RegisterFunctionUtilityInterpolant> > ScatteringMapUtiltiyEnum;
    //
    typedef std::map<std::string, Enum> MapStringToEnum;

    // Scattering* CreateScattering(const std::string&, Particle&, Utility&);

    Scattering* CreateScattering(const std::string&, Particle&, Utility&, const InterpolationDef&);
    Scattering* CreateScattering(const Enum, Particle&, Utility&, const InterpolationDef&);

    Scattering* CreateScattering(const std::string&, Particle&, Utility&);
    Scattering* CreateScattering(const Enum, Particle&, Utility&);

    // Scattering* CreateScattering(const Enum, Particle&, Utility&);

    Enum GetEnumFromString(const std::string&);

    static ScatteringFactory& Get()
    {
        static ScatteringFactory instance;
        return instance;
    }

    private:
    ScatteringFactory();
    ~ScatteringFactory();

    void Register(const std::string& name, Enum);

    std::vector<Enum> registerd_enum;
    std::vector<std::string> registerd_str;
    MapStringToEnum map_string_to_enum;
    // void Register(const std::string& name, Enum, RegisterFunction);
    // void RegisterUtility(const std::string& name, Enum, std::pair<RegisterFunctionUtility, RegisterFunctionUtilityInterpolant> RegisterFunctionUtility);

    // std::map<std::string, Scattering* (*)(void)> scattering_map;
    // ScatteringMapString scattering_map_str_;
    // ScatteringMapEnum scattering_map_enum_;
    //
    // ScatteringMapUtiltiyString scattering_map_utility_str_;
    // ScatteringMapUtiltiyEnum scattering_map_utility_enum_;
    //
};

} /*  PROPOSAL */

