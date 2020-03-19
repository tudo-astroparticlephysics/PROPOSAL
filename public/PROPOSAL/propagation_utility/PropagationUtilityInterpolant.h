
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

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include <functional>
#include <memory>
#include <string>

using std::string;

namespace PROPOSAL {
class Integral;
class Interpolant;
struct InterpolationDef;

class UtilityInterpolant : public UtilityIntegral {
    Interpolant1DBuilder builder1d;
    Interpolant1DBuilder builder_diff;
    std::unique_ptr<Interpolant> interpolant_;
    std::unique_ptr<Interpolant> interpolant_diff_;
    double stored_result_;
    double low_;

    void BuildTables(const string, size_t);
public:
    UtilityInterpolant(std::function<double(double)>);

    virtual double Calculate(double ei, double ef, double rnd);
    virtual double GetUpperLimit(double ei, double rnd);

    static InterpolationDef utility_interpolation_def;
};
} // namespace PROPOSAL
