
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

#include <map>
#include <string>
#include <memory>
#include <sstream>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

    class CrossSection;
    class Compton;
    struct ParticleDef;
    class EnergyCutSettings;
    class Medium;

    class ComptonFactory
    {
    public:
        // --------------------------------------------------------------------- //
        // Declare usable enums
        // --------------------------------------------------------------------- //

        enum Enum
        {
            Fail = 0,
            None,
            KleinNishina
        };

        struct Definition
        {
            Definition()
                    : parametrization(None)
                    , multiplier(1.0)
            {
            }

            bool operator==(const ComptonFactory::Definition& def) const
            {
                if (parametrization != def.parametrization)
                    return false;
                else if (multiplier != def.multiplier)
                    return false;

                return true;
            }

            bool operator!=(const ComptonFactory::Definition& def) const
            {
                return !(*this == def);
            }

            friend std::ostream& operator<<(std::ostream& os, PROPOSAL::ComptonFactory::Definition const& definition)
            {
                std::stringstream ss;
                ss << " Compton Definition (" << &definition << ") ";
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
                Compton*(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier)>
                RegisterFunction;

        typedef std::map<std::string, RegisterFunction> ComptonMapString;
        typedef std::map<Enum, RegisterFunction> ComptonMapEnum;
        typedef Helper::Bimap<std::string, Enum> BimapStringEnum;

        // --------------------------------------------------------------------- //
        // Most general creation
        // --------------------------------------------------------------------- //

        CrossSection* CreateCompton(const ParticleDef&,
                                           std::shared_ptr<const Medium>,
                                           const EnergyCutSettings&,
                                           const Definition&) const;

        CrossSection* CreateCompton(const ParticleDef&,
                                           std::shared_ptr<const Medium>,
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

        static ComptonFactory& Get()
        {
            static ComptonFactory instance;
            return instance;
        }

    private:
        ComptonFactory();
        ~ComptonFactory();

        // ----------------------------------------------------------------------------
        /// @brief Register Compton parametrizations
        ///
        /// @param name
        /// @param Enum
        /// @param RegisterFunction
        // ----------------------------------------------------------------------------
        void Register(const std::string& name, Enum, RegisterFunction);

        ComptonMapString compton_map_str_;
        ComptonMapEnum compton_map_enum_;
        BimapStringEnum string_enum_;
    };

} // namespace PROPOSAL
