
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

namespace PROPOSAL {
namespace crosssection {
    class Annihilation : public Parametrization<Component> {
    public:
        Annihilation() = default;
        ~Annihilation() override = default;

        double GetLowerEnergyLim(ParticleDef const&) const noexcept final;

        KinematicLimits GetKinematicLimits(
            ParticleDef const&, Component const&, double) const final;
    };

    template <> struct ParametrizationName<Annihilation> {
        static constexpr auto value = "annihilation";
    };

    struct AnnihilationHeitler : public Annihilation {
        AnnihilationHeitler() = default;

        std::unique_ptr<Parametrization<Component>> clone() const final;

        double DifferentialCrossSection(
            ParticleDef const&, Component const&, double, double) const final;
    };

    template <> struct ParametrizationName<AnnihilationHeitler> {
        static constexpr auto value = "Heitler";
    };

    template <> struct ParametrizationId<AnnihilationHeitler> {
        static constexpr size_t value = 1000000012;
    };

} // namespace crosssection
} // namespace PROPOSAL
