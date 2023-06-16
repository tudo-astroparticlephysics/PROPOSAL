
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
#include "PROPOSAL/math/Vector3D.h"
#include <memory>
#include <utility>

namespace PROPOSAL {
namespace multiple_scattering {
    struct ScatteringOffset {
        ScatteringOffset() : sx(0.), sy(0.), tx(0.), ty(0.) {};
        friend std::ostream& operator<<(std::ostream& os, const ScatteringOffset& offset) {
            os << "sx: " << offset.sx << ", sy: " << offset.sy << ", tx: " << offset.tx << ", ty: " << offset.ty  << '\n';
            return os;
        };
        double sx;
        double sy;
        double tx;
        double ty;
    };

    std::pair<Cartesian3D, Cartesian3D> ScatterInitialDirection(
            const Vector3D&, const ScatteringOffset&);

    class Parametrization {
    protected:
        double mass;

    public:
        Parametrization() = default;
        Parametrization(double mass);
        virtual ~Parametrization() = default;

        virtual std::unique_ptr<Parametrization> clone() const = 0;

        virtual bool compare(const Parametrization&) const = 0;
        virtual void print(std::ostream&) const = 0;

        virtual ScatteringOffset CalculateRandomAngle(double grammage,
            double ei, double ef, const std::array<double, 4>& rnd)
            = 0;

        // Calculate a single scattering angle
        virtual double CalculateScatteringAngle(double grammage, double ei, double ef, double rnd) = 0;

        // Calculate a scattering angle in space, i.e. sqrt(theta_1**2 + theta_2**), where theta_1 and theta_2 are
        // independent scattering angles
        virtual double CalculateScatteringAngle2D(double grammage, double ei, double ef, double rnd1, double rnd2) = 0;

        virtual bool operator==(const Parametrization& scattering) const;
        friend std::ostream& operator<<(std::ostream&, Parametrization const&);
    };
} // namespace multiple_scattering
} // namespace PROPOSAL
