
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <boost/bimap.hpp>
#include <functional>

#include <map>
#include <string>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

class CrossSection;
class RealPhoton;
class ShadowEffect;
class Photonuclear;
struct ParticleDef;
class EnergyCutSettings;
class Medium;

class PhotonuclearFactory
{
public:
    // --------------------------------------------------------------------- //
    // Declare usable enums
    // --------------------------------------------------------------------- //

    enum Enum
    {
        None = 0,
        Zeus,
        BezrukovBugaev,
        Rhode,
        Kokoulin,
        AbramowiczLevinLevyMaor91,
        AbramowiczLevinLevyMaor97,
        ButkevichMikhailov,
        RenoSarcevicSu
    };

    enum Shadow
    {
        ShadowNone,
        ShadowDuttaRenoSarcevicSeckel,
        ShadowButkevichMikhailov
    };

    struct Definition
    {
        Definition()
            : parametrization(AbramowiczLevinLevyMaor97)
            , shadow(ShadowButkevichMikhailov)
            , hard_component(true)
            , multiplier(1.0)
        {
        }

        bool operator==(const PhotonuclearFactory::Definition& def) const
        {
            if (parametrization != def.parametrization)
                return false;
            else if (shadow != def.shadow)
                return false;
            else if (hard_component != def.hard_component)
                return false;
            else if (multiplier != def.multiplier)
                return false;

            return true;
        }

        bool operator!=(const PhotonuclearFactory::Definition& def) const
        {
            return !(*this == def);
        }

        Enum parametrization;
        Shadow shadow;
        bool hard_component;
        double multiplier;
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef std::function<ShadowEffect*(void)> RegisterShadowEffectFunction;

    typedef std::function<Photonuclear*(const ParticleDef&,
                                          const Medium&,
                                          const EnergyCutSettings&,
                                          double multiplier,
                                          bool hard_component)>
        RegisterRealPhotonFunction;

    typedef std::function<Photonuclear*(const ParticleDef&,
                                          const Medium&,
                                          const EnergyCutSettings&,
                                          double multiplier,
                                          const ShadowEffect&)>
        RegisterQ2Function;

    typedef std::function<Photonuclear*(const ParticleDef&,
                                          const Medium&,
                                          const EnergyCutSettings&,
                                          double multiplier,
                                          const ShadowEffect&,
                                          InterpolationDef)>
        RegisterQ2FunctionInterpolant;

    typedef std::map<std::string, RegisterShadowEffectFunction> PhotoShadowEffectMapString;
    typedef std::map<Shadow, RegisterShadowEffectFunction> PhotoShadowEffectMapEnum;

    typedef std::map<std::string, RegisterRealPhotonFunction> PhotoRealPhotonMapString;
    typedef std::map<Enum, RegisterRealPhotonFunction> PhotoRealPhotonMapEnum;

    typedef std::map<std::string, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> > PhotoQ2MapString;
    typedef std::map<Enum, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> > PhotoQ2MapEnum;

    typedef boost::bimap<std::string, Enum> BimapStringEnum;
    typedef boost::bimap<std::string, Shadow> BimapStringShadowEnum;

    // --------------------------------------------------------------------- //
    // Shadow effect
    // --------------------------------------------------------------------- //

    ShadowEffect* CreateShadowEffect(const std::string&);
    ShadowEffect* CreateShadowEffect(const Shadow&);

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreatePhotonuclear(const ParticleDef&,
                                     const Medium&,
                                     const EnergyCutSettings&,
                                     const Definition&) const;

    CrossSection* CreatePhotonuclear(const ParticleDef&,
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

    // ----------------------------------------------------------------------------
    /// @brief Enum string conversation for shadow effect
    // ----------------------------------------------------------------------------
    Shadow GetShadowEnumFromString(const std::string&);

    // ----------------------------------------------------------------------------
    /// @brief enum to string conversation for shadow effect
    // ----------------------------------------------------------------------------
    std::string GetStringFromShadowEnum(const Shadow&);

    // --------------------------------------------------------------------- //
    // Singleton pattern
    // --------------------------------------------------------------------- //

    static PhotonuclearFactory& Get()
    {
        static PhotonuclearFactory instance;
        return instance;
    }

private:
    PhotonuclearFactory();
    ~PhotonuclearFactory();

    // ----------------------------------------------------------------------------
    /// @brief Register ShadowEffect used for Q2 integration parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param std::pair
    // ----------------------------------------------------------------------------
    void RegisterShadowEffect(const std::string& name, const Shadow&, RegisterShadowEffectFunction);

    // ----------------------------------------------------------------------------
    /// @brief Register Photonuclear parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterRealPhotonFunction
    // ----------------------------------------------------------------------------
    void RegisterRealPhoton(const std::string& name, Enum, RegisterRealPhotonFunction);

    // ----------------------------------------------------------------------------
    /// @brief Register Photonuclear parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterQ2Function
    // ----------------------------------------------------------------------------
    void RegisterQ2(const std::string& name, Enum, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant>);

    PhotoShadowEffectMapString photo_shadow_map_str_;
    PhotoShadowEffectMapEnum photo_shadow_map_enum_;

    PhotoRealPhotonMapString photo_real_map_str_;
    PhotoRealPhotonMapEnum photo_real_map_enum_;

    PhotoQ2MapString photo_q2_map_str_;
    PhotoQ2MapEnum photo_q2_map_enum_;

    BimapStringEnum string_enum_;
    BimapStringShadowEnum string_shadow_enum_;
};

} // namespace PROPOSAL
