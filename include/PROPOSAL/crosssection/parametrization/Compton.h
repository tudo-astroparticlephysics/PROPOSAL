
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

    struct ComptonKleinNishina : public Compton {
        ComptonKleinNishina() = default;

        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double energy, double v) const final;
    };
} // namespace crosssection
} // namespace PROPOSAL

/* namespace detail { */
/*     template <typename T1, typename T2> */
/*     auto define_dndx_integral( */
/*         crosssection::ComptonKleinNishina param, T1 p_def, T2 target) */
/*     { */
/*         return [param, p_def, target](Integral& i, double E, double v_min, */
/*                    double v_max, double rate) { */
/*             auto t_min = std::log(1. - v_max); */
/*             auto t_max = std::log(1. - v_min); */
/*             auto dNdx = [&param, &p_def, &target, E](double t) { */
/*                 return exp(t) */
/*                     * param.DifferentialCrossSection( */
/*                         p_def, target, E, 1 - exp(t)); */
/*             }; */
/*             return i.IntegrateWithRandomRatio(t_min, t_max, dNdx, 3, rate);
 */
/*         }; */
/*     } */

/*     template <typename T1, typename T2> */
/*     auto define_dndx_upper_lim( */
/*         crosssection::ComptonKleinNishina param, T1 p_def, T2 target) */
/*     { */
/*         return [param, p_def, target](Integral& i, double E, double v_min, */
/*                    double v_max, double rnd) { */
/*             auto t_min = std::log(1. - v_min); */
/*             auto t_max = std::log(1. - v_max); */
/*             auto dNdx = [&param, &p_def, &target, E](double t) { */
/*                 return exp(t) */
/*                     * param.DifferentialCrossSection( */
/*                         p_def, target, E, 1 - exp(t)); */
/*             }; */
/*             i.IntegrateWithRandomRatio(t_min, t_max, dNdx, 3, rnd); */
/*             return 1. - std::exp(i.GetUpperLimit()); */
/*         }; */
/*     } */

/*     inline auto define_dedx_integral(crosssection::ComptonKleinNishina param,
 */
/*         ParticleDef const& p_def, Component const& comp, */
/*         EnergyCutSettings const& cut) */
/*     { */
/*         return [param, p_def, comp, cut](Integral& i, double E) { */
/*             auto physical_lim = param.GetKinematicLimits(p_def, comp, E); */
/*             auto v_cut = cut.GetCut(physical_lim, E); */
/*             auto t_min = std::log(1. - v_cut); */
/*             auto t_max = std::log(1. */
/*                 -
 * std::get<crosssection::Parametrization::V_MIN>(physical_lim)); */
/*             auto dEdx = [&param, &p_def, &comp, E](double t) { */
/*                 return exp(t) */
/*                     * param.FunctionToDEdxIntegral(p_def, comp, E, 1 -
 * exp(t)); */
/*             }; */
/*             return i.Integrate(t_min, t_max, dEdx, 2); */
/*         }; */
/*     } */

/*     inline auto define_de2dx_integral(crosssection::ComptonKleinNishina
 * param, */
/*         ParticleDef const& p_def, Component const& comp, */
/*         EnergyCutSettings const& cut) */
/*     { */
/*         return [param, p_def, comp, cut](Integral& i, double E) { */
/*             auto physical_lim = param.GetKinematicLimits(p_def, comp, E); */
/*             auto v_cut = cut.GetCut(physical_lim, E); */
/*             auto t_min = std::log(1. - v_cut); */
/*             auto t_max = std::log(1. */
/*                 -
 * std::get<crosssection::Parametrization::V_MIN>(physical_lim)); */
/*             auto dEdx = [&param, &p_def, &comp, E](double t) { */
/*                 return exp(t) */
/*                     * param.FunctionToDE2dxIntegral(p_def, comp, E, 1 -
 * exp(t)); */
/*             }; */
/*             return i.Integrate(t_min, t_max, dEdx, 2); */
/*         }; */
/*     } */
}

// Factory pattern functions

/* namespace crosssection { */
/*     template <typename P, typename M> */
/*     using compt_func_ptr = cross_t_ptr<P, M> (*)( */
/*         P, M, std::shared_ptr<const EnergyCutSettings>, bool); */

/*     template <typename Param, typename P, typename M> */
/*     cross_t_ptr<P, M> create_compt(P p_def, M medium, */
/*         std::shared_ptr<const EnergyCutSettings> cuts, bool interpol) */
/*     { */
/*         auto param = Param(); */
/*         return make_crosssection(param, p_def, medium, cuts, interpol); */
/*     } */

/*     template <typename P, typename M> */
/*     static std::map<std::string, compt_func_ptr<P, M>> compt_map = { */
/*         { "kleinnishina", create_compt<ComptonKleinNishina, P, M> }, */
/*     }; */

/*     template <typename P, typename M> */
/*     cross_t_ptr<P, M> make_compton(P p_def, M medium, */
/*         std::shared_ptr<const EnergyCutSettings> cuts, bool interpol, */
/*         const std::string& param_name) */
/*     { */
/*         std::string name = param_name; */
/*         std::transform( */
/*             param_name.begin(), param_name.end(), name.begin(), ::tolower);
 */
/*         auto it = compt_map<P, M>.find(name); */
/*         if (it == compt_map<P, M>.end()) */
/*             throw std::logic_error("Unknown parametrization for compton"); */

/*         return it->second(p_def, medium, cuts, interpol); */
/*     } */

/*     template <typename P, typename M> */
/*     cross_t_ptr<P, M> make_compton(P p_def, M medium, */
/*         std::shared_ptr<const EnergyCutSettings> cuts, bool interpol, */
/*         const nlohmann::json& config) */
/*     { */
/*         if (!config.contains("parametrization")) */
/*             throw std::logic_error("No parametrization passed for compton");
 */
/*         std::string param_name = config["parametrization"]; */

/*         return make_compton(p_def, medium, cuts, interpol, param_name); */
/*     } */
/* } // namespace crosssection */
