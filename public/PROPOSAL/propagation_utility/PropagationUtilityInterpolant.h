
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

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

namespace PROPOSAL {
class Integral;
class UtilityIntegral;
} // namespace PROPOSAL

#define UTILITY_INTERPOLANT_DEC(cls)                                                                                   \
    class UtilityInterpolant##cls : public UtilityInterpolant                                                          \
    {                                                                                                                  \
    public:                                                                                                            \
        UtilityInterpolant##cls(const Utility&, InterpolationDef);                                                     \
                                                                                                                       \
        UtilityInterpolant##cls(const Utility&, const UtilityInterpolant##cls&);                                       \
        UtilityInterpolant##cls(const UtilityInterpolant##cls&);                                                       \
        UtilityDecorator* clone(const Utility& utility) const { return new UtilityInterpolant##cls(utility, *this); }  \
                                                                                                                       \
        ~UtilityInterpolant##cls();                                                                                    \
                                                                                                                       \
        double Calculate(double ei, double ef, double rnd);                                                            \
        double GetUpperLimit(double ei, double rnd);                                                                   \
                                                                                                                       \
    private:                                                                                                           \
        UtilityInterpolant##cls& operator=(const UtilityInterpolant##cls&);                                            \
                                                                                                                       \
        double BuildInterpolant(double, UtilityIntegral&, Integral&);                                                  \
        void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);                   \
    };

namespace PROPOSAL {

class Interpolant;

class UtilityInterpolant : public UtilityDecorator
{
public:
    UtilityInterpolant(const Utility&, InterpolationDef);

    // Copy constructors
    UtilityInterpolant(const Utility&, const UtilityInterpolant&);
    UtilityInterpolant(const UtilityInterpolant&);
    virtual UtilityDecorator* clone(const Utility&) const = 0;

    virtual ~UtilityInterpolant();

    // Methods
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd);

protected:
    UtilityInterpolant& operator=(const UtilityInterpolant&); // Undefined & not allowed

    virtual bool compare(const UtilityDecorator&) const;

    virtual double BuildInterpolant(double, UtilityIntegral&, Integral&)                                = 0;
    virtual void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points) = 0;

    double stored_result_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    InterpolationDef interpolation_def_;
};

class UtilityInterpolantInteraction : public UtilityInterpolant
{
public:
    UtilityInterpolantInteraction(const Utility&, InterpolationDef);

    // Copy constructors
    UtilityInterpolantInteraction(const Utility&, const UtilityInterpolantInteraction&);
    UtilityInterpolantInteraction(const UtilityInterpolantInteraction&);
    virtual UtilityInterpolant* clone(const Utility& utility) const
    {
        return new UtilityInterpolantInteraction(utility, *this);
    }

    virtual ~UtilityInterpolantInteraction();

    // Methods
    double Calculate(double ei, double ef, double rnd);
    double GetUpperLimit(double ei, double rnd);

private:
    UtilityInterpolantInteraction& operator=(const UtilityInterpolantInteraction&); // Undefined & not allowed

    virtual bool compare(const UtilityDecorator&) const;

    double BuildInterpolant(double, UtilityIntegral&, Integral&);
    void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);

    double big_low_;
    double up_;
};

class UtilityInterpolantDecay : public UtilityInterpolant
{
public:
    UtilityInterpolantDecay(const Utility&, InterpolationDef);

    // Copy constructors
    UtilityInterpolantDecay(const Utility&, const UtilityInterpolantDecay&);
    UtilityInterpolantDecay(const UtilityInterpolantDecay&);
    virtual UtilityInterpolant* clone(const Utility& utility) const
    {
        return new UtilityInterpolantDecay(utility, *this);
    }

    virtual ~UtilityInterpolantDecay();

    // Methods
    double Calculate(double ei, double ef, double rnd);
    double GetUpperLimit(double ei, double rnd);

private:
    UtilityInterpolantDecay& operator=(const UtilityInterpolantDecay&); // Undefined & not allowed

    virtual bool compare(const UtilityDecorator&) const;

    double BuildInterpolant(double, UtilityIntegral&, Integral&);
    void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);

    double big_low_;
    double up_;
};

UTILITY_INTERPOLANT_DEC(Displacement)
UTILITY_INTERPOLANT_DEC(Time)
UTILITY_INTERPOLANT_DEC(ContRand)
UTILITY_INTERPOLANT_DEC(Scattering)

} // namespace PROPOSAL
