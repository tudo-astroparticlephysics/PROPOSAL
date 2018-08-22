
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

#include <map>
#include <string>
#include <vector>

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class Medium;
class Utility;
struct InterpolationDef;

class ScatteringFactory
{
public:
    enum Enum
    {
        None = 0,
        HighlandIntegral,
        Moliere,
        Highland,
        NoScattering
    };

    typedef boost::bimap<std::string, Enum> BimapStringEnum;

    Scattering* CreateScattering(const std::string&, Particle&, const Utility&, const InterpolationDef&);
    Scattering* CreateScattering(const Enum, Particle&, const Utility&, const InterpolationDef&);

    Scattering* CreateScattering(const std::string&, Particle&, const Utility&);
    Scattering* CreateScattering(const Enum, Particle&, const Utility&);

    // ----------------------------------------------------------------------------
    /// @brief string to enum conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);

    // ----------------------------------------------------------------------------
    /// @brief enum to string conversation for photo parametrizations
    // ----------------------------------------------------------------------------
    std::string GetStringFromEnum(const Enum&);

    static ScatteringFactory& Get()
    {
        static ScatteringFactory instance;
        return instance;
    }

private:
    ScatteringFactory();
    ~ScatteringFactory();

    void Register(const std::string& name, const Enum);

    std::vector<Enum> registerd_enum;
    std::vector<std::string> registerd_str;
    BimapStringEnum string_enum_;
};

} // namespace PROPOSAL
