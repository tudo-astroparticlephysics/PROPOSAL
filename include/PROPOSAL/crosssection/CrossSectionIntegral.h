
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

/* #include "PROPOSAL/EnergyCutSettings.h" */
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/Integral.h"
/* #include "PROPOSAL/medium/Medium.h" */
/* #include "PROPOSAL/particle/ParticleDef.h" */

namespace PROPOSAL {
template <typename Param, typename P, typename M>
class CrossSectionIntegral : public crosssection_t<P, M> {
    using param_t = typename std::decay<Param>::type;

public:
    CrossSectionIntegral(Param _param, P _p_def, M _medium,
        std::shared_ptr<const EnergyCutSettings> _cut)
        : crosssection_t<P, M>(_param, _p_def, _medium,
            detail::build_dndx(
                typename param_t::base_param_t::component_wise {}, false,
                _medium, _param, _p_def, _cut),
            detail::build_dedx(
                typename param_t::base_param_t::component_wise {}, false,
                _param, _p_def, _medium, _cut),
            detail::build_de2dx(
                typename param_t::base_param_t::component_wise {}, false,
                _param, _p_def, _medium, _cut))
    {
        if (typename param_t::only_stochastic {} == true and _cut != nullptr) {
            throw std::invalid_argument(
                "CrossSections of parametrizations that are only stochastic do "
                "not use"
                "EnergyCuts. Pass a nullptr as an EnergyCut instead.");
        }
    }

    inline double CalculatedNdx(double energy,
        std::shared_ptr<const Component> comp_ptr = nullptr) override
    {
        return this->CalculatedNdx_impl(energy,
            typename param_t::base_param_t::component_wise {}, comp_ptr);
    }
    inline double CalculateStochasticLoss(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate) override
    {
        assert(rate > 0);
        return this->CalculateStochasticLoss_impl(comp, energy, rate,
            typename param_t::base_param_t::component_wise {},
            typename param_t::base_param_t::only_stochastic {});
    }
};
} // namespace PROPOSAL
