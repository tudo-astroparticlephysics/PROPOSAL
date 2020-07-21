
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

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

using std::add_lvalue_reference;
using std::decay;

namespace PROPOSAL {

using Components::Component;

using rates_t = std::unordered_map<const Component*, double>;

template <class P, class M> struct CrossSection {
    CrossSection() = default;
    virtual ~CrossSection() = default;

    virtual double CalculatedEdx(double) = 0;
    virtual double CalculatedE2dx(double) = 0;
    virtual rates_t CalculatedNdx(double) = 0;
    virtual double CalculateStochasticLoss(const Component&, double, double)
        = 0;

    virtual size_t GetHash() const noexcept = 0;
    virtual double GetLowerEnergyLim() const = 0;
    virtual InteractionType GetInteractionType() const noexcept = 0;
};

template <typename P, typename M>
using crosssection_t
    = CrossSection<typename std::decay<P>::type, typename std::decay<M>::type>;

template <typename P, typename M>
using crosssection_list_t = std::vector<std::shared_ptr<crosssection_t<P, M>>>;

} // namespace PROPOSAL
