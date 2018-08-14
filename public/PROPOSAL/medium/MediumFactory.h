
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

    typedef std::function<Medium*(double)> RegisterFunction;
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
