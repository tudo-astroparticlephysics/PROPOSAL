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
#include "PROPOSAL/medium/density_distr/density_distr.h"

namespace PROPOSAL {

class Density_homogeneous : public Density_distr {
   public:
    Density_homogeneous();
    Density_homogeneous(const Density_homogeneous&);
    Density_homogeneous(double correction_factor);

    // ~Density_homogeneous(){};

    bool compare(const Density_distr& dens_distr) const override;

    Density_distr* clone() const override {
        return new Density_homogeneous(*this);
    };

    double Correct(const Vector3D& xi,
                   const Vector3D& direction,
                   double res,
                   double distance_to_border) const override;
    double Integrate(const Vector3D& xi,
                     const Vector3D& direction,
                     double res) const override;
    double Calculate(const Vector3D& xi,
                     const Vector3D& direction,
                     double distance) const override;
    double Evaluate(const Vector3D& xi) const override;

    double GetCorrectionfactor() const { return correction_factor_; }

   private:
    double correction_factor_;
};
}  // namespace PROPOSAL
