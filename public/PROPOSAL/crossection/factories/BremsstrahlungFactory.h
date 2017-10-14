
#pragma once

#include <boost/function.hpp>

#include <map>
#include <string>

#include "PROPOSAL/methods.h"

namespace  PROPOSAL
{

class CrossSection;
class Bremsstrahlung;
struct ParticleDef;
class EnergyCutSettings;

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

    typedef std::map<std::string, RegisterFunction > BremsstrahlungMapString;
    typedef std::map<Enum, RegisterFunction > BremsstrahlungMapEnum;
    typedef std::map<std::string, Enum> MapStringToEnum;

    // --------------------------------------------------------------------- //
    // Create functions
    // --------------------------------------------------------------------- //

    CrossSection* CreateBremsstrahlungIntegral(const std::string&,
                                       const ParticleDef& particle_def,
                                       const Medium& medium,
                                       const EnergyCutSettings& cuts,
                                       double multiplier,
                                       bool lpm) const;

    CrossSection* CreateBremsstrahlungIntegral(const Enum,
                                       const ParticleDef& particle_def,
                                       const Medium& medium,
                                       const EnergyCutSettings& cuts,
                                       double multiplier,
                                       bool lpm) const;

    CrossSection* CreateBremsstrahlungInterpolant(const std::string&,
                                       const ParticleDef& particle_def,
                                       const Medium& medium,
                                       const EnergyCutSettings& cuts,
                                       double multiplier,
                                       bool lpm,
                                       InterpolationDef) const;

    CrossSection* CreateBremsstrahlungInterpolant(const Enum,
                                       const ParticleDef& particle_def,
                                       const Medium& medium,
                                       const EnergyCutSettings& cuts,
                                       double multiplier,
                                       bool lpm,
                                       InterpolationDef) const;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateBremsstrahlung(const Enum,
                                       const ParticleDef&,
                                       const Medium&,
                                       const EnergyCutSettings&,
                                       double multiplier,
                                       bool lpm,
                                       bool interpolate,
                                       InterpolationDef = InterpolationDef()) const;

    CrossSection* CreateBremsstrahlung(const ParticleDef&,
                                       const Medium&,
                                       const EnergyCutSettings&,
                                       const Definition&,
                                       bool interpolate,
                                       InterpolationDef = InterpolationDef()) const;

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

