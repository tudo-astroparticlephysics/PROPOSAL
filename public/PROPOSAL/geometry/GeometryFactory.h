
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

#include "PROPOSAL/geometry/Geometry.h"

namespace PROPOSAL {

class GeometryFactory
{
public:
    enum Enum
    {
        Sphere = 0,
        Box,
        Cylinder
    };

    struct Definition
    {
        Definition()
            : shape(Sphere)
            , position()
            , inner_radius(0.0)
            , radius(0.0)
            , width(0.0)
            , height(0.0)
            , depth(0.0)
        {
        }

        Enum shape;
        Vector3D position;
        double inner_radius;
        double radius;
        double width;
        double height;
        double depth;
    };

    typedef std::function<Geometry*(void)> RegisterFunction;
    typedef std::map<std::string, RegisterFunction> GeometryMapString;
    typedef std::map<Enum, RegisterFunction> GeometryMapEnum;

    void Register(const std::string&, const Enum&, RegisterFunction);

    Geometry* CreateGeometry(const std::string&);
    Geometry* CreateGeometry(const Enum&);
    Geometry* CreateGeometry(const Definition&);

    static GeometryFactory& Get()
    {
        static GeometryFactory instance;
        return instance;
    }

private:
    GeometryFactory();
    ~GeometryFactory();

    GeometryMapString geometry_map_str;
    GeometryMapEnum geometry_map_enum;
};

} // namespace PROPOSAL
