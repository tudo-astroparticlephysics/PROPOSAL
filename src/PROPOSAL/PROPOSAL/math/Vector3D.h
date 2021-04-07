
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
#include <array>
#include <ostream>
namespace PROPOSAL {
class Cartesian3D;
class Spherical3D;
class Vector3D {
public:
    Vector3D() : coordinates({0, 0, 0}) {};
    Vector3D(std::array<double, 3> val) : coordinates(val) {};
    virtual ~Vector3D() = default;

    virtual bool operator==(const Vector3D&) const;
    virtual bool operator!=(const Vector3D&) const;
    double& operator[](size_t idx);
    const double& operator[](size_t idx) const;
    friend std::ostream& operator<<(std::ostream&, const Vector3D&);

    virtual double magnitude() const = 0;
    virtual void normalize() = 0;

    void SetCoordinates(std::array<double, 3> val) {coordinates = val;}
    void SetCoordinates(double, double, double);
    virtual std::array<double, 3> GetSphericalCoordinates() const = 0;
    virtual std::array<double, 3> GetCartesianCoordinates() const = 0;

protected:
    std::array<double, 3> coordinates;
    virtual void print(std::ostream&) const = 0;

};
} // namespace PROPOSAL
