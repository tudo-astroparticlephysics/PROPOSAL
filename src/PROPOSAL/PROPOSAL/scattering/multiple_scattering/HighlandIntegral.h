
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

#include "PROPOSAL/scattering/multiple_scattering/Highland.h"

namespace PROPOSAL {
    class Displacement;
    struct CrossSectionBase;
    class UtilityIntegral;
}

namespace PROPOSAL {
namespace multiple_scattering {
    class HighlandIntegral : public Highland {
        std::unique_ptr<UtilityIntegral> highland_integral;

        double CalculateTheta0(double, double, double) final;
        inline double Integral(Displacement&, double);

    public:
        HighlandIntegral(const ParticleDef& p, Medium const& m,
            std::shared_ptr<Displacement> disp, std::false_type);

        HighlandIntegral(const ParticleDef& p, Medium const& m,
            std::shared_ptr<Displacement> disp, std::true_type);
    };
} // namespace multiple_scattering

std::unique_ptr<multiple_scattering::Parametrization> make_highland_integral(
        ParticleDef const& p, Medium const& m,
        std::shared_ptr<Displacement> disp, bool interpol = false);

    std::unique_ptr<multiple_scattering::Parametrization> make_highland_integral(
            ParticleDef const& p, Medium const& m,
            std::vector<std::shared_ptr<CrossSectionBase>> const& c,
            bool interpol = false);
} // namespace PROPOSAL
