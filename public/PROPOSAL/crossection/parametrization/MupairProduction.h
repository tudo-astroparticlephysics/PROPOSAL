
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

#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/Integral.h"

#define MUPAIR_PARAM_INTEGRAL_DEC(param)                                       \
    class Mupair##param : public MupairProductionRhoIntegral {                 \
    public:                                                                    \
        Mupair##param(const ParticleDef&, const component_list&);              \
                                                                               \
        double FunctionToIntegral(double energy, double v, double r);          \
    };

namespace PROPOSAL {

class MupairProduction : public Parametrization {
protected:
    Integral drho_integral_;

public:
    MupairProduction(const ParticleDef&, const component_list&);
    virtual ~MupairProduction() = default;

    virtual double DifferentialCrossSection(double energy, double v) = 0;
    virtual double FunctionToIntegral(double energy, double v, double rho) = 0;
    double Calculaterho(double energy, double v, double rnd1, double rnd2);

    KinematicLimits GetKinematicLimits(double energy);
};

class MupairProductionRhoIntegral : public MupairProduction {
    Integral integral_;
public:
    MupairProductionRhoIntegral(const ParticleDef&, const component_list&);
    virtual ~MupairProductionRhoIntegral() = default;

    virtual double DifferentialCrossSection(double energy, double v);
};

MUPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)

#undef MUPAIR_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
