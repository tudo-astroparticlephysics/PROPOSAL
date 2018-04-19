
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

#define UTILITY_INTEGRAL_DEC(cls)                                                                                      \
    class UtilityIntegral##cls : public UtilityIntegral                                                                \
    {                                                                                                                  \
    public:                                                                                                            \
        UtilityIntegral##cls(const Utility&);                                                                          \
                                                                                                                       \
        UtilityIntegral##cls(const Utility&, const UtilityIntegral##cls&);                                             \
        UtilityIntegral##cls(const UtilityIntegral##cls&);                                                             \
        virtual ~UtilityIntegral##cls();                                                                               \
                                                                                                                       \
        virtual UtilityDecorator* clone(const Utility& utility) const                                                  \
        {                                                                                                              \
            return new UtilityIntegral##cls(utility, *this);                                                           \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToIntegral(double energy);                                                                      \
        virtual double Calculate(double ei, double ef, double rnd);                                                    \
                                                                                                                       \
    private:                                                                                                           \
        UtilityDecorator& operator=(const UtilityDecorator&);                                                          \
    };

namespace PROPOSAL {

class UtilityIntegral : public UtilityDecorator
{
public:
    UtilityIntegral(const Utility&);

    // Copy constructors
    UtilityIntegral(const Utility&, const UtilityIntegral&);
    UtilityIntegral(const UtilityIntegral&);
    virtual UtilityDecorator* clone(const Utility&) const = 0;

    virtual ~UtilityIntegral();

    // Methods
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd);

protected:
    UtilityIntegral& operator=(const UtilityIntegral&); // Undefined & not allowed

    virtual bool compare(const UtilityDecorator&) const;

    Integral integral_;
};

class UtilityIntegralDisplacement : public UtilityIntegral
{
public:
    UtilityIntegralDisplacement(const Utility&);

    // Copy constructors
    UtilityIntegralDisplacement(const Utility&, const UtilityIntegralDisplacement&);
    UtilityIntegralDisplacement(const UtilityIntegralDisplacement&);
    virtual UtilityIntegral* clone(const Utility& utility) const
    {
        return new UtilityIntegralDisplacement(utility, *this);
    }

    virtual ~UtilityIntegralDisplacement();

    // Methods
    virtual double Calculate(double ei, double ef, double rnd);

private:
    UtilityDecorator& operator=(const UtilityDecorator&); // Undefined & not allowed
};

UTILITY_INTEGRAL_DEC(Interaction)
UTILITY_INTEGRAL_DEC(Decay)
UTILITY_INTEGRAL_DEC(Time)
UTILITY_INTEGRAL_DEC(ContRand)
UTILITY_INTEGRAL_DEC(Scattering)

#undef UTILITY_INTEGRAL_DEC

} // namespace PROPOSAL
