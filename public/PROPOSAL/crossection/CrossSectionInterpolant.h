
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

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {
class CrossSectionInterpolant : public CrossSectionIntegral {
    size_t hash_interpol_def;

protected:
    unique_ptr<Interpolant> dedx_interpolant_;
    unique_ptr<Interpolant> de2dx_interpolant_;
    vector<unique_ptr<Interpolant>> dndx_interpolants_;

    double logarithm_trafo(double, double, double) const;

    virtual vector<unique_ptr<Interpolant>> init_dndx_interpolation(
        const InterpolationDef&);
    unique_ptr<Interpolant> init_dedx_interpolation(const InterpolationDef&);
    unique_ptr<Interpolant> init_de2dx_interpolation(const InterpolationDef&);

public:
    CrossSectionInterpolant(unique_ptr<Parametrization>&&,
        shared_ptr<const EnergyCutSettings>, const InterpolationDef&);

    double CalculatedEdx(double) override;
    double CalculatedE2dx(double) override;
    double CalculatedNdx(double, double = 1) override;

    size_t GetHash() const override;
};
} // namespace PROPOSAL
