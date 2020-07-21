
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

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

namespace PROPOSAL {
static std::map<const Geometry_Type, std::shared_ptr<Geometry>> Geometry_Map{
    { Geometry_Type::SPHERE,
        std::make_shared<Sphere>(Vector3D(0, 0, 0), 0.f, 0.f) },
    { Geometry_Type::BOX,
        std::make_shared<Box>(Vector3D(0, 0, 0), 0.f, 0.f, 0.f) },
    { Geometry_Type::CYLINDER,
        std::make_shared<Cylinder>(Vector3D(0, 0, 0), 0.f, 0.f, 0.f) },
};
} // namespace PROPOSAL

namespace PROPOSAL {
std::shared_ptr<Geometry> CreateGeometry(std::string name);
std::shared_ptr<Geometry> CreateGeometry(const nlohmann::json& config);
std::shared_ptr<Geometry> CreateGeometry(Geometry_Type type);
} // namespace PROPOSAL
