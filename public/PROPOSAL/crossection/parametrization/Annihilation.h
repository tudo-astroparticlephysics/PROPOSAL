
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

#include <cmath>
#include <fstream>
#include <functional>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace crosssection {
class Annihilation : public Parametrization {
public:
    Annihilation();
    virtual ~Annihilation() = default;

    using component_wise = std::true_type;
    using only_stochastic = std::true_type;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    tuple<double, double> GetKinematicLimits(const ParticleDef&, const Component&, double) const noexcept override;
};

struct AnnihilationHeitler : public Annihilation {
    AnnihilationHeitler();

    using base_param_t = Annihilation;

    double DifferentialCrossSection(
        const ParticleDef&, const Component&, double, double) override;
};
} // namespace crosssection
} // namespace PROPOSAL
