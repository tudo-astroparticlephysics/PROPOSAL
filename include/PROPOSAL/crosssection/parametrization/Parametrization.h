
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

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include <tuple>

using PROPOSAL::Components::Component;

namespace PROPOSAL {
class Medium;
enum class InteractionType;
namespace crosssection {
    struct Parametrization {
        const InteractionType interaction_type;
        const std::string name;
        size_t hash;

        Parametrization(InteractionType, const std::string&);
        virtual ~Parametrization() = default;

        virtual double DifferentialCrossSection(
            const ParticleDef&, const Component&, double, double) const = 0;

        enum { V_MIN, V_MAX };
        virtual std::tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Component&, double) const = 0;

        inline double FunctionToDEdxIntegral(ParticleDef const& p_def,
            const Component& comp, double energy, double v) const
        {
            return v * DifferentialCrossSection(p_def, comp, energy, v);
        }

        inline double FunctionToDE2dxIntegral(ParticleDef const& p_def,
            const Component& comp, double energy, double v) const
        {
            return v * v * DifferentialCrossSection(p_def, comp, energy, v);
        }

        virtual double GetLowerEnergyLim(ParticleDef const&) const noexcept = 0;

        inline size_t GetHash() const noexcept { return hash; };
    };
} // namespace crosssection
} // namespace PROPOSAL
