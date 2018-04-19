
#pragma once

#include <boost/bimap.hpp>
#include <boost/function.hpp>

#include <map>

namespace PROPOSAL {

class Medium;

class MediumFactory
{
public:
    enum Enum
    {
        Water = 0,
        Ice,
        Salt,
        StandardRock,
        FrejusRock,
        Iron,
        Hydrogen,
        Lead,
        Copper,
        Uranium,
        Air,
        Paraffin,
        AntaresWater
    };

    struct Definition
    {
        Definition()
            : type(Water)
            , density_correction(1.0)
        {
        }

        Enum type;
        double density_correction;
    };

    typedef boost::function<Medium*(double)> RegisterFunction;
    typedef std::map<std::string, RegisterFunction> MediumMapString;
    typedef std::map<Enum, RegisterFunction> MediumMapEnum;
    typedef boost::bimap<std::string, Enum> BimapStringEnum;

    void Register(const std::string& name, const Enum&, RegisterFunction);
    // void Register(const Enum&, RegisterFunction);

    Medium* CreateMedium(const std::string&, double density_correction = 1.0);
    Medium* CreateMedium(const Enum&, double density_correction = 1.0);
    Medium* CreateMedium(Definition);

    // ----------------------------------------------------------------------------
    /// @brief string to enum conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);

    // ----------------------------------------------------------------------------
    /// @brief enum to string conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    std::string GetStringFromEnum(const Enum&);

    static MediumFactory& Get()
    {
        static MediumFactory instance;
        return instance;
    }

private:
    MediumFactory();
    ~MediumFactory();

    MediumMapString medium_map_str;
    MediumMapEnum medium_map_enum;

    BimapStringEnum string_enum_;
};

} // namespace PROPOSAL
