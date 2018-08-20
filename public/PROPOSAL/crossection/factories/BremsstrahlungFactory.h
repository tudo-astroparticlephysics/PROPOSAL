
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
        None = 0,
        PetrukhinShestakov,
        KelnerKokoulinPetrukhin,
        CompleteScreening,
        AndreevBezrukovBugaev,
        SandrockSoedingreksoRhode
    };

    struct Definition
    {
        Definition()
            : parametrization(KelnerKokoulinPetrukhin)
            , lpm_effect(true)
            , multiplier(1.0)
        {
        }

        bool operator==(const BremsstrahlungFactory::Definition& def) const
        {
            if (parametrization != def.parametrization)
                return false;
            else if (lpm_effect != def.lpm_effect)
                return false;
            else if (multiplier != def.multiplier)
                return false;

            return true;
        }

        bool operator!=(const BremsstrahlungFactory::Definition& def) const
        {
            return !(*this == def);
        }

        Enum parametrization;
        bool lpm_effect;
        double multiplier;
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef std::function<
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
