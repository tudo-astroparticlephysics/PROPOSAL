
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
} // namespace PROPOSAL;

namespace PROPOSAL {
namespace crosssection {

    class WeakInteraction : public Parametrization {
    public:
        WeakInteraction() = default;
        virtual ~WeakInteraction() = default;

        double GetLowerEnergyLim(ParticleDef const&) const final;
        KinematicLimits GetKinematicLimits(
            ParticleDef const&, Component const&, double) const final;
    };

    template <> struct ParametrizationName<WeakInteraction> {
        static constexpr char value[36] = "weak_interaction";
    };

    template <> struct ParametrizationId<WeakInteraction> {
        static constexpr size_t value = 1000000009;
    };

    class WeakCooperSarkarMertsch : public WeakInteraction {
        using Interpolant_t = std::shared_ptr<Interpolant>;
        std::pair<Interpolant_t, Interpolant_t> interpolants_particle;
        std::pair<Interpolant_t, Interpolant_t> interpolants_antiparticle;

    public:
        WeakCooperSarkarMertsch();

        double DifferentialCrossSection(
            ParticleDef const&, Component const&, double, double) const final;
    };

    template <> struct ParametrizationName<WeakCooperSarkarMertsch> {
        static constexpr char value[36] = "weak_cooper_sarkar_mertsch";
    };
} // namespace crosssection
} // namespace PROPOSAL

/* // Factory pattern functions */

/* template <typename P, typename M> */
/* using weak_func_ptr = cross_t_ptr<P, M> (*)(P, M, bool); */

/* template <typename Param, typename P, typename M> */
/* cross_t_ptr<P, M> create_weak(P p_def, M medium, bool interpol) */
/* { */
/*     auto param = Param(); */
/*     return make_crosssection(param, p_def, medium, nullptr, interpol); */
/* } */

/* template <typename P, typename M> */
/* static std::map<std::string, weak_func_ptr<P, M>> weak_map = { */
/*     { "coopersarkarmertsch", create_weak<WeakCooperSarkarMertsch, P, M> }
 */
/* }; */

/* template <typename P, typename M> */
/* cross_t_ptr<P, M> make_weakinteraction( */
/*     P p_def, M medium, bool interpol, const std::string& param_name) */
/* { */
/*     std::string name = param_name; */
/*     std::transform( */
/*         param_name.begin(), param_name.end(), name.begin(), ::tolower);
 */
/*     auto it = weak_map<P, M>.find(name); */
/*     if (it == weak_map<P, M>.end()) */
/*         throw std::logic_error( */
/*             "Unknown parametrization for weak interaction"); */

/*     return it->second(p_def, medium, interpol); */
/* } */

/* template <typename P, typename M> */
/* cross_t_ptr<P, M> make_weakinteraction( */
/*     P p_def, M medium, bool interpol, const nlohmann::json& config) */
/* { */
/*     if (!config.contains("parametrization")) */
/*         throw std::logic_error( */
/*             "No parametrization passed for weak interaction"); */
/*     std::string param_name = config["parametrization"]; */

/*     return make_weakinteraction(p_def, medium, interpol, param_name); */
/* } */
