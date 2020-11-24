
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
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include <functional>
#include <memory>
#include <string>

namespace PROPOSAL {
class Integral;
class Interpolant;
struct InterpolationDef;

class UtilityInterpolant : public UtilityIntegral {
    std::unique_ptr<Interpolant> interpolant_;
    std::pair<double, double> upper_limit;

    // maybe interpolate function to integral will give a performance boost.
    // in general this function should be underfrequently called
    // Interpolant1DBuilder builder_diff;
    // std::unique_ptr<Interpolant> interpolant_diff_;

public:
    UtilityInterpolant(std::function<double(double)>, double);
    void BuildTables(const std::string, size_t, Interpolant1DBuilder::Definition, bool = false);

    double Calculate(double, double);
    double GetUpperLimit(double, double);
};
} // namespace PROPOSAL
