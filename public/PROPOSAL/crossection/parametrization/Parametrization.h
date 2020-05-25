
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
#include "PROPOSAL/medium/Components.h"
#include <iostream>
#include <type_traits>

using PROPOSAL::Components::Component;
using std::string;

namespace PROPOSAL {
class ParticleDef;
class Medium;
enum class InteractionType;

class Parametrization {
public:
    const InteractionType interaction_type;
    const string name;

    Parametrization(InteractionType, const string&);
    virtual ~Parametrization() = default;

    struct KinematicLimits {
        double vMin;
        double vMax;

        KinematicLimits(double, double);
    };

    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double, double)
        = 0;
    virtual double FunctionToDNdxIntegral(
        const ParticleDef&, const Component&, double, double);
    virtual double FunctionToDEdxIntegral(
        const ParticleDef&, const Component&, double, double);
    virtual double FunctionToDE2dxIntegral(
        const ParticleDef&, const Component&, double, double);
    virtual KinematicLimits GetKinematicLimits(
        const ParticleDef&, const Component&, double)
        = 0;

    virtual double GetLowerEnergyLim(const ParticleDef&) const { return 0; };
    virtual size_t GetHash() const;
};
} // namespace PROPOSAL

namespace PROPOSAL {
template <typename Param>
double integrate_dndx(Integral& integral, Param&& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max)
{
    auto dNdx = [&param, &p_def, &comp, energy](double v) {
        return param.FunctionToDNdxIntegral(p_def, comp, energy, v);
    };
    return integral.Integrate(v_min, v_max, dNdx, 4);
}

template <typename Param>
double calculate_upper_lim_dndx(Integral& integral, Param&& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max, double rnd)
{
    auto dNdx = [&param, &p_def, &comp, energy](double v) {
        return param.FunctionToDNdxIntegral(p_def, comp, energy, v);
    };
    integral.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rnd);
    return integral.GetUpperLimit();
}

template <typename Param>
double integrate_dedx(Integral& integral, Param&& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max)
{
    auto dEdx = [&param, &p_def, &comp, energy](double v) {
        return param.FunctionToDEdxIntegral(p_def, comp, energy, v);
    };
    return integral.Integrate(v_min, v_max, dEdx, 2);
}

template <typename Param>
double integrate_de2dx(Integral& integral, Param&& param,
    const ParticleDef& p_def, const Component& comp, double energy,
    double v_min, double v_max)
{
    auto dE2dx = [&param, &p_def, &comp, energy](double v) {
        return param.FunctionToDE2dxIntegral(p_def, comp, energy, v);
    };
    return integral.Integrate(v_min, v_max, dE2dx, 2);
}
} // namespace PROPOSAL
