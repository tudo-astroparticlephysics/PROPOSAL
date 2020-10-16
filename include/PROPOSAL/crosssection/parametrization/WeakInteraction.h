
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
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/math/Interpolant.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace PROPOSAL {
class Interpolant;

namespace crosssection {
    class WeakInteraction : public Parametrization {
    public:
        WeakInteraction();
        ~WeakInteraction() = default;

        using only_stochastic = std::true_type;
        using component_wise = std::true_type;

        double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
        std::tuple<double, double> GetKinematicLimits(const ParticleDef&,
            const Component&, double) const noexcept override;
    };

    class WeakCooperSarkarMertsch : public WeakInteraction {
        std::unordered_map<bool,
            std::tuple<std::shared_ptr<Interpolant>,
                std::shared_ptr<const Interpolant>>>
            interpolant_;
        std::tuple<Interpolant, Interpolant> BuildContribution(
            bool is_decayable) const;

    public:
        WeakCooperSarkarMertsch();

        using base_param_t = WeakInteraction;

        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double, double) const override;
    };

// Factory pattern functions

template <typename P, typename M>
using weak_func_ptr = cross_t_ptr<P, M>(*)(P, M, bool);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_weak(P p_def, M medium, bool interpol) {
    auto param = Param();
    return make_crosssection(param, p_def, medium, nullptr, interpol);
}

template<typename P, typename M>
static std::map<std::string, weak_func_ptr<P, M>> weak_map = {
        {"CooperSarkarMertsch", create_weak<WeakCooperSarkarMertsch, P, M>}
};

template<typename P, typename M>
cross_t_ptr<P, M> make_weakinteraction(P p_def, M medium, bool interpol,
                                    const std::string& param_name){

    auto it = weak_map<P, M>.find(param_name);
    if (it == weak_map<P, M>.end())
        throw std::logic_error("Unknown parametrization for weak interaction");

    return it->second(p_def, medium, interpol);
}

template<typename P, typename M>
cross_t_ptr<P, M> make_weakinteraction(P p_def, M medium, bool interpol,
                                    const nlohmann::json& config){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for weak interaction");
    std::string param_name = config["parametrization"];

    return make_weakinteraction(p_def, medium, interpol, param_name);
}

} // namespace crosssection
} // namespace PROPOSAL
