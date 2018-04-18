
#pragma once

// #include <boost/function.hpp>
#include <boost/bimap.hpp>

#include <map>
#include <string>
#include <vector>

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class Medium;
class Utility;
struct InterpolationDef;

class ScatteringFactory
{
public:
    enum Enum
    {
        HighlandIntegral = 0,
        Moliere,
        Highland,
        NoScattering
    };

    typedef boost::bimap<std::string, Enum> BimapStringEnum;

    Scattering* CreateScattering(const std::string&, Particle&, const Utility&, const InterpolationDef&);
    Scattering* CreateScattering(const Enum, Particle&, const Utility&, const InterpolationDef&);

    Scattering* CreateScattering(const std::string&, Particle&, const Utility&);
    Scattering* CreateScattering(const Enum, Particle&, const Utility&);

    // ----------------------------------------------------------------------------
    /// @brief string to enum conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);

    // ----------------------------------------------------------------------------
    /// @brief enum to string conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    std::string GetStringFromEnum(const Enum&);

    static ScatteringFactory& Get()
    {
        static ScatteringFactory instance;
        return instance;
    }

private:
    ScatteringFactory();
    ~ScatteringFactory();

    void Register(const std::string& name, const Enum);

    std::vector<Enum> registerd_enum;
    std::vector<std::string> registerd_str;
    BimapStringEnum string_enum_;
};

} // namespace PROPOSAL
