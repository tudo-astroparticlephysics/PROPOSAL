
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

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using PROPOSAL::Components::Component;

namespace PROPOSAL {
class Compton : public Parametrization {
public:
    Compton();
    virtual ~Compton() = default;

    using only_stochastic = std::false_type;
    using component_wise = std::true_type;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double) const noexcept override;
};

struct ComptonKleinNishina : public Compton {
    ComptonKleinNishina() = default;
    using base_param_t = Compton;

    double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v) const;
};

template <>
double integrate_dndx(Integral&, Compton&, const ParticleDef&, const Component&,
    double, double, double, double);

template <>
double calculate_upper_lim_dndx(Integral&, Compton&, const ParticleDef&,
    const Component&, double, double, double, double);

template <>
double integrate_dedx(Integral&, Compton&, const ParticleDef&, const Component&,
    double, double, double);

template <>
double integrate_de2dx(Integral&, Compton&, const ParticleDef&,
    const Component&, double, double, double);
} // namespace PROPOSAL
