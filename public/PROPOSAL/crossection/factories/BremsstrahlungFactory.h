
#pragma once

#include <boost/bimap.hpp>
#include <boost/function.hpp>

#include <map>
#include <string>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

class CrossSection;
class Bremsstrahlung;
struct ParticleDef;
class EnergyCutSettings;
class Medium;

class BremsstrahlungFactory
{
public:
    // --------------------------------------------------------------------- //
    // Declare usable enums
    // --------------------------------------------------------------------- //

    enum Enum
    {
        PetrukhinShestakov = 0,
        KelnerKokoulinPetrukhin,
        CompleteScreening,
        AndreevBezrukovBugaev
    };

    struct Definition
    {
        Definition()
            : parametrization(KelnerKokoulinPetrukhin)
            , lpm_effect(true)
            , multiplier(1.0)
        {
        }

        Enum parametrization;
        bool lpm_effect;
        double multiplier;
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef boost::function<
        Bremsstrahlung*(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, bool lpm)>
        RegisterFunction;

    typedef std::map<std::string, RegisterFunction> BremsstrahlungMapString;
    typedef std::map<Enum, RegisterFunction> BremsstrahlungMapEnum;
    typedef boost::bimap<std::string, Enum> BimapStringEnum;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateBremsstrahlung(const ParticleDef&,
                                       const Medium&,
                                       const EnergyCutSettings&,
                                       const Definition&) const;

    CrossSection* CreateBremsstrahlung(const ParticleDef&,
                                       const Medium&,
                                       const EnergyCutSettings&,
                                       const Definition&,
                                       InterpolationDef) const;

    // ----------------------------------------------------------------------------
    /// @brief string to enum conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);

    // ----------------------------------------------------------------------------
    /// @brief enum to string conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    std::string GetStringFromEnum(const Enum&);

    // --------------------------------------------------------------------- //
    // Singleton pattern
    // --------------------------------------------------------------------- //

    static BremsstrahlungFactory& Get()
    {
        static BremsstrahlungFactory instance;
        return instance;
    }

private:
    BremsstrahlungFactory();
    ~BremsstrahlungFactory();

    // ----------------------------------------------------------------------------
    /// @brief Register Bremsstrahlung parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterFunction
    // ----------------------------------------------------------------------------
    void Register(const std::string& name, Enum, RegisterFunction);

    BremsstrahlungMapString bremsstrahlung_map_str_;
    BremsstrahlungMapEnum bremsstrahlung_map_enum_;
    BimapStringEnum string_enum_;
};

} // namespace PROPOSAL
