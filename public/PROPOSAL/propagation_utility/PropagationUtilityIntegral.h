
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
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

namespace PROPOSAL {
class UtilityIntegral : public UtilityDecorator {
public:
    UtilityIntegral(CrossSectionList, const ParticleDef&);

    virtual double Calculate(double ei, double ef, double rnd);
    virtual double GetUpperLimit(double ei, double rnd);
    virtual double FunctionToIntegral(double energy) = 0;

protected:
    Integral integral_;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityIntegralDisplacement : public UtilityIntegral {
public:
    UtilityIntegralDisplacement(CrossSectionList, const ParticleDef&);
    double FunctionToIntegral(double energy) override;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityIntegralInteraction : public UtilityIntegral {
public:
    UtilityIntegralInteraction(CrossSectionList, const ParticleDef&);
    double FunctionToIntegral(double energy) override;

protected:
    UtilityIntegralDisplacement displacement_;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityIntegralDecay : public UtilityIntegral {
public:
    UtilityIntegralDecay(CrossSectionList, const ParticleDef&);
    double FunctionToIntegral(double energy) override;
    double Calculate(double ei, double ef, double rnd) override;

protected:
    UtilityIntegralDisplacement displacement_;
    double lifetime;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityIntegralTime : public UtilityIntegral {
public:
    UtilityIntegralTime(CrossSectionList, const ParticleDef&);
    double FunctionToIntegral(double energy) override;

protected:
    UtilityIntegralDisplacement displacement_;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityIntegralContRand : public UtilityIntegral {
public:
    UtilityIntegralContRand(CrossSectionList, const ParticleDef&);
    double FunctionToIntegral(double energy) override;

protected:
    UtilityIntegralDisplacement displacement_;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityIntegralScattering : public UtilityIntegral {
public:
    UtilityIntegralScattering(CrossSectionList, const ParticleDef&);
    double FunctionToIntegral(double energy) override;

protected:
    UtilityIntegralDisplacement displacement_;
};
} // namespace PROPOSAL
