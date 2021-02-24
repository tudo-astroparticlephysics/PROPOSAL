
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
    class Compton : public Parametrization<Component> {
    public:
        Compton() = default;
        virtual ~Compton() = default;

        double GetLowerEnergyLim(ParticleDef const&) const noexcept final;
        KinematicLimits GetKinematicLimits(
            ParticleDef const&, Component const&, double) const final;
    };

    template <> struct ParametrizationName<Compton> {
        static constexpr char value[36] = "compton";
    };

    struct ComptonKleinNishina : public Compton {
        ComptonKleinNishina() = default;

        std::unique_ptr<Parametrization<Component>> clone() const final;

        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double energy, double v) const final;
    };

    template <> struct ParametrizationName<ComptonKleinNishina> {
        static constexpr auto value = "KleinNishina";
    };

    template <> struct ParametrizationId<ComptonKleinNishina> {
        static constexpr size_t value = 1000000010;
    };

} // namespace crosssection
} // namespace PROPOSAL
