
#pragma once

#include <boost/function.hpp>

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace  PROPOSAL
{

class CrossSection;

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

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef boost::function<
        Parametrization*(const ParticleDef&, const Medium&, const EnergyCutSettings&, Parametrization::Definition)>
        RegisterFunction;

    typedef std::map<std::string, RegisterFunction > BremsstrahlungMapString;
    typedef std::map<Enum, RegisterFunction > BremsstrahlungMapEnum;
    typedef std::map<std::string, Enum> MapStringToEnum;

    // --------------------------------------------------------------------- //
    // Create functions
    // --------------------------------------------------------------------- //

    Parametrization* CreateParametrization(const std::string&,
                                           const ParticleDef& particle_def,
                                           const Medium& medium,
                                           const EnergyCutSettings& cuts,
                                           Parametrization::Definition def) const;

    Parametrization* CreateParametrization(const Enum,
                                           const ParticleDef& particle_def,
                                           const Medium& medium,
                                           const EnergyCutSettings& cuts,
                                           Parametrization::Definition def) const;

    CrossSection* CreateBremsstrahlung(const std::string&,
                                       const ParticleDef& particle_def,
                                       const Medium& medium,
                                       const EnergyCutSettings& cuts,
                                       Parametrization::Definition def,
                                       bool interpolate = true) const;

    CrossSection* CreateBremsstrahlung(const Enum,
                                       const ParticleDef& particle_def,
                                       const Medium& medium,
                                       const EnergyCutSettings& cuts,
                                       Parametrization::Definition def,
                                       bool interpolate = true) const;


    // ----------------------------------------------------------------------------
    /// @brief Enum string conversation
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);


    // ----------------------------------------------------------------------------
    /// @brief Register Bremsstrahlung parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterFunction
    // ----------------------------------------------------------------------------
    void Register(const std::string& name, Enum, RegisterFunction);

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


    BremsstrahlungMapString bremsstrahlung_map_str_;
    BremsstrahlungMapEnum bremsstrahlung_map_enum_;
    MapStringToEnum map_string_to_enum;
};

} /*  PROPOSAL */

