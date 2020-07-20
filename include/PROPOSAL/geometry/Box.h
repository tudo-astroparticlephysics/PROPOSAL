
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

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/json.hpp"


namespace PROPOSAL {

class Box : public Geometry
{
public:
    Box();
    Box(const Vector3D position, double x, double y, double z);
    Box(const Box&);
    Box(const nlohmann::json& config);

    std::shared_ptr<const Geometry> create() const override { return std::shared_ptr<const Geometry>( new Box(*this) ); };
    void swap(Geometry&) override;

    virtual ~Box() {}

    // Operators
    Box& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;

    // Getter & Setter
    double GetX() const { return x_; }
    double GetY() const { return y_; }
    double GetZ() const { return z_; }

    void SetX(double x) { x_ = x; };
    void SetY(double y) { y_ = y; };
    void SetZ(double z) { z_ = z; };

private:
    bool compare(const Geometry&) const override;
    void print(std::ostream&) const override;

    double x_; //!< width of box in x-direction
    double y_; //!< width of box in y-direction
    double z_; //!< width of box in z-direction
};

} // namespace PROPOSAL
