
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
#include <iostream>
#include <vector>
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/Spline.h"
#include "PROPOSAL/medium/density_distr/density_polynomial.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namespace PROPOSAL {
class Density_splines : public Density_distr {
   public:
    Density_splines(const Axis&, const Spline&);
    Density_splines(const Density_splines&);
    Density_splines() {
        delete axis_;
        delete spline_;
        delete integrated_spline_;
    }

    bool compare(const Density_distr& dens_distr) const override;

    Density_distr* clone() const override {
        return new Density_splines(*this);
    };

    double Correct(const Vector3D& xi,
                   const Vector3D& direction,
                   double res,
                   double distance_to_border) const override;
    double Integrate(const Vector3D& xi, const Vector3D& direction, double l) const override;
    double Evaluate(const Vector3D& xi) const override;
    double Calculate(const Vector3D& xi,
                     const Vector3D& direction,
                     double distance) const override;

    double Helper_function(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double l) const;
    double helper_function(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double l) const;

   protected:
    Spline* spline_;
    Spline* integrated_spline_;

    std::function<double(double)> density_distribution;
    std::function<double(double)> antiderived_density_distribution;
};
}  // namespace PROPOSAL
