
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
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <iostream>
#include <tuple>
#include <type_traits>

using PROPOSAL::Components::Component;

namespace PROPOSAL {
class Medium;
enum class InteractionType;
namespace crosssection {

    struct Parametrization {
        const InteractionType interaction_type;
        const std::string name;

        Parametrization(InteractionType, const std::string&);
        virtual ~Parametrization() = default;

        virtual double DifferentialCrossSection(
            const ParticleDef&, const Component&, double, double) const = 0;

        enum { V_MIN, V_MAX };
        virtual std::tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Component&, double) const = 0;

        inline double FunctionToDEdxIntegral(const ParticleDef& p_def,
            const Component& comp, double energy, double v) const
        {
            return v * DifferentialCrossSection(p_def, comp, energy, v);
        }
        inline double FunctionToDE2dxIntegral(const ParticleDef& p_def,
            const Component& comp, double energy, double v) const
        {
            return v * v * DifferentialCrossSection(p_def, comp, energy, v);
        }
        virtual double GetLowerEnergyLim(const ParticleDef&) const noexcept = 0;
        virtual size_t GetHash() const noexcept;
    };
} // namespace crosssection
} // namespace PROPOSAL

namespace PROPOSAL {
template <typename P, typename M>
double integrate_dndx(Integral& integral, P&& param, const ParticleDef& p_def,
    const M& medium, double energy, double v_min, double v_max, double rnd = 0)
{
    auto dNdx = [&param, &p_def, &medium, energy](double v) {
        return param.DifferentialCrossSection(p_def, medium, energy, v);
    };
    return integral.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rnd);
}

template <typename P, typename M>
double calculate_upper_lim_dndx(Integral& integral, P&& param,
    const ParticleDef& p_def, M&& medium, double energy, double v_min,
    double v_max, double rnd)
{
    auto dNdx = [&param, &p_def, &medium, energy](double v) {
        return param.DifferentialCrossSection(p_def, medium, energy, v);
    };
    return integral.IntegrateWithRandomRatio(v_min, v_max, dNdx, 4, rnd);
    return integral.GetUpperLimit();
}

template <typename Param, typename M>
double integrate_dedx(Integral& integral, Param&& param,
    const ParticleDef& p_def, M& medium, double energy, double v_min,
    double v_max)
{
    auto dEdx = [&param, &p_def, &medium, energy](double v) {
        return param.FunctionToDEdxIntegral(p_def, medium, energy, v);
    };
    return integral.Integrate(v_min, v_max, dEdx, 2);
}

template <typename Param, typename M>
double integrate_de2dx(Integral& integral, Param&& param,
    const ParticleDef& p_def, M& medium, double energy, double v_min,
    double v_max)
{
    auto dE2dx = [&param, &p_def, &medium, energy](double v) {
        return param.FunctionToDE2dxIntegral(p_def, medium, energy, v);
    };
    return integral.Integrate(v_min, v_max, dE2dx, 2);
}

} // namespace PROPOSAL
