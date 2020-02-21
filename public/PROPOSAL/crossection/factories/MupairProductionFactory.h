
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
#include <sstream>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

class CrossSection;
class MupairProduction;

struct ParticleDef;
class Medium;
class EnergyCutSettings;

class MupairProductionFactory
{
public:
    enum Enum
    {
        Fail = 0,
        None,
        KelnerKokoulinPetrukhin,
    };

    struct Definition
    {
        Definition()
            : parametrization(None)
            , multiplier(1.0)
            , particle_output(true)
        {
        }

        bool operator==(const MupairProductionFactory::Definition& def) const
        {
            if (parametrization != def.parametrization)
                return false;
            else if (multiplier != def.multiplier)
                return false;
            else if (particle_output != def.particle_output)
                return false;

            return true;
        }

        bool operator!=(const MupairProductionFactory::Definition& def) const
        {
            return !(*this == def);
        }

        friend std::ostream& operator<<(std::ostream& os, PROPOSAL::MupairProductionFactory::Definition const& definition)
        {
            std::stringstream ss;
            ss << " MuPair Production Definition (" << &definition << ") ";
            os << Helper::Centered(60, ss.str()) << '\n';

            os << "Parametrization: " << definition.parametrization << std::endl;
            os << "Multiplier: " << definition.multiplier << std::endl;
            os << "Particle Output: " << definition.particle_output << std::endl;

            os << Helper::Centered(60, "");
            return os;
        }

        Enum parametrization;
        double multiplier;
        bool particle_output;
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef std::function<MupairProduction*(const ParticleDef&,
                                          std::shared_ptr<const Medium>,
                                          const EnergyCutSettings&,
                                          double multiplier,
                                          bool particle_output)>
        RegisterFunction;

    typedef std::function<MupairProduction*(const ParticleDef&,
                                          std::shared_ptr<const Medium>,
                                          const EnergyCutSettings&,
                                          double multiplier,
                                          bool particle_output,
                                          InterpolationDef)>
        RegisterFunctionInterpolant;

    typedef std::map<std::string, std::pair<RegisterFunction, RegisterFunctionInterpolant> > MupairMapString;
    typedef std::map<Enum, std::pair<RegisterFunction, RegisterFunctionInterpolant> > MupairMapEnum;

    typedef Helper::Bimap<std::string, Enum> BimapStringEnum;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateMupairProduction(const ParticleDef&,
                                        std::shared_ptr<const Medium>,
                                        const EnergyCutSettings&,
                                        const Definition&) const;

    CrossSection* CreateMupairProduction(const ParticleDef&,
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

    static MupairProductionFactory& Get()
    {
        static MupairProductionFactory instance;
        return instance;
    }

private:
    MupairProductionFactory();
    ~MupairProductionFactory();

    // ----------------------------------------------------------------------------
    /// @brief Register MuPairProduction parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterFunction
    // ----------------------------------------------------------------------------
    void Register(const std::string& name, Enum, std::pair<RegisterFunction, RegisterFunctionInterpolant>);

    MupairMapString mupair_map_str_;
    MupairMapEnum mupair_map_enum_;
    BimapStringEnum string_enum_;
};

} // namespace PROPOSAL
