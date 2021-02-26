
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

#include "PROPOSAL/math/Integral.h"
#include <string>

namespace PROPOSAL {
class UtilityIntegral {
    Integral integral;

protected:
    double lower_lim;
    std::function<double(double)> FunctionToIntegral;
    size_t hash;

public:
    UtilityIntegral() = default;
    UtilityIntegral(std::function<double(double)>, double, size_t);
    virtual ~UtilityIntegral() = default;

    virtual void BuildTables(const std::string, size_t, bool) {};

    virtual double Calculate(double, double);
    virtual double GetUpperLimit(double, double);

    virtual size_t GetHash() const { return hash; }
};
} // namespace PROPOSAL
