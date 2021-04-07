
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

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <memory>
#include <utility>

namespace PROPOSAL {
class Interpolant;
class Component;
} // namespace PROPOSAL;

namespace PROPOSAL {
namespace crosssection {

    struct WeakInteraction : public Parametrization<Component> {
    public:
        WeakInteraction() = default;
        virtual ~WeakInteraction() = default;

        double GetLowerEnergyLim(ParticleDef const&) const noexcept final;
        KinematicLimits GetKinematicLimits(
            ParticleDef const&, Component const&, double) const final;
    };

    template <> struct ParametrizationName<WeakInteraction> {
        static constexpr char value[36] = "weak_interaction";
    };

    template <> struct ParametrizationId<WeakInteraction> {
        static constexpr size_t value = 1000000009;
    };

    struct WeakCooperSarkarMertsch : public WeakInteraction {
        using Interpolant_t = std::shared_ptr<Interpolant>;
        std::pair<Interpolant_t, Interpolant_t> interpolants_particle;
        std::pair<Interpolant_t, Interpolant_t> interpolants_antiparticle;

    public:
        WeakCooperSarkarMertsch();

        std::unique_ptr<Parametrization<Component>> clone() const final;

        double DifferentialCrossSection(
            ParticleDef const&, Component const&, double, double) const final;
    };

    template <> struct ParametrizationName<WeakCooperSarkarMertsch> {
        static constexpr auto value = "CooperSarkarMertsch";
    };

    template <> struct ParametrizationId<WeakCooperSarkarMertsch> {
        static constexpr size_t value = 1000000009;
    };
} // namespace crosssection
} // namespace PROPOSAL
