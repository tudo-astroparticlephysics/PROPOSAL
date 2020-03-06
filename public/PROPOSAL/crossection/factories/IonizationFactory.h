
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

#include <functional>
#include <memory>

#include <map>
#include <string>
#include <iostream>
#include <sstream>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

class CrossSection;
class Ionization;
struct ParticleDef;
class Medium;
class EnergyCutSettings;

class IonizationFactory
{
public:
    enum Enum
    {
        Fail = 0,
        None,
        BetheBlochRossi,
        IonizBergerSeltzerBhabha,
        IonizBergerSeltzerMoller,
    };

    struct Definition
    {
        Definition()
            : parametrization(BetheBlochRossi)
            , multiplier(1.0)
        {
        }

        bool operator==(const IonizationFactory::Definition& def) const
        {
            if(parametrization != def.parametrization)
                return false;
            else if(multiplier != def.multiplier)
                return false;
            return true;
        }

        bool operator!=(const IonizationFactory::Definition& def) const
        {
            return !(*this == def);
        }

        friend std::ostream& operator<<(std::ostream& os, PROPOSAL::IonizationFactory::Definition const& definition)
        {
            std::stringstream ss;
            ss << " Ionization Definition (" << &definition << ") ";
            os << Helper::Centered(60, ss.str()) << '\n';

            os << "Parametrization: " << definition.parametrization << std::endl;
            os << "Multiplier: " << definition.multiplier << std::endl;

            os << Helper::Centered(60, "");
            return os;
        }

        Enum parametrization;
        double multiplier;
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef std::function<
            Ionization*(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier)>
            RegisterFunction;

    typedef std::map<std::string, RegisterFunction> IonizMapString;
    typedef std::map<Enum, RegisterFunction> IonizMapEnum;
        typedef Helper::Bimap<std::string, Enum> BimapStringEnum;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateIonization(const ParticleDef&,
                                   std::shared_ptr<const Medium>,
                                   const EnergyCutSettings&,
                                   const Definition&) const;

    CrossSection* CreateIonization(const ParticleDef&,
                                   std::shared_ptr<const Medium>,
                                   const EnergyCutSettings&,
                                   const Definition&,
                                   InterpolationDef) const;

    // ----------------------------------------------------------------------------
    /// @brief string to enum conversation for ioniz parametrizations
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);

    // ----------------------------------------------------------------------------
    /// @brief enum to string conversation for ioniz parametrizations
    // ----------------------------------------------------------------------------
    std::string GetStringFromEnum(const Enum&);

    // --------------------------------------------------------------------- //
    // Singleton pattern
    // --------------------------------------------------------------------- //

    static IonizationFactory& Get()
    {
        static IonizationFactory instance;
        return instance;
    }

private:
    IonizationFactory();
    ~IonizationFactory();

    // ----------------------------------------------------------------------------
    /// @brief Register Ionization parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterFunction
    // ----------------------------------------------------------------------------
    void Register(const std::string& name, Enum, RegisterFunction);

    IonizMapString ioniz_map_str_;
    IonizMapEnum ioniz_map_enum_;
    BimapStringEnum string_enum_;
};

} // namespace PROPOSAL
