
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
#include <functional>

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/math/Integral.h"

#define MUPAIR_PARAM_INTEGRAL_DEC(param)                                       \
    struct Mupair##param : public MupairProductionRhoIntegral {                \
        Mupair##param();                                                       \
        using base_param_t = MupairProduction;                                 \
                                                                               \
        double FunctionToIntegral(const ParticleDef&, const Component&,        \
            double energy, double v, double r) const;                          \
    };

namespace PROPOSAL {
namespace crosssection {

    class MupairProduction : public Parametrization {

    protected:
        Integral drho_integral_;

    public:
        MupairProduction();
        virtual ~MupairProduction() = default;
        using only_stochastic = std::false_type;
        using component_wise = std::true_type;

        virtual double FunctionToIntegral(const ParticleDef&, const Component&,
            double energy, double v, double r) const = 0;

        double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
        std::tuple<double, double> GetKinematicLimits(const ParticleDef&,
            const Component&, double) const noexcept override;
    };

    class MupairProductionRhoIntegral : public MupairProduction {

    public:
        MupairProductionRhoIntegral();
        ~MupairProductionRhoIntegral() = default;

        virtual double DifferentialCrossSection(
            const ParticleDef&, const Component&, double, double) const;
    };

    MUPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)
// Factory pattern functions

template <typename P, typename M>
using mupair_func_ptr = cross_t_ptr<P, M>(*)(P, M, std::shared_ptr<const
        EnergyCutSettings>, bool);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_mupair(P p_def, M medium,std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol) {
    auto param = Param();
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template<typename P, typename M>
static std::map<std::string, mupair_func_ptr<P, M>> mupair_map = {
        {"KelnerKokoulinPetrukhin", create_mupair<MupairKelnerKokoulinPetrukhin, P, M>},
};

template<typename P, typename M>
cross_t_ptr<P, M> make_mupairproduction(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, const std::string& param_name){

    auto it = mupair_map<P, M>.find(param_name);
    if (it == mupair_map<P, M>.end())
        throw std::logic_error("Unknown parametrization for mupairproduction");

    return it->second(p_def, medium, cuts, interpol);
}

template<typename P, typename M>
cross_t_ptr<P, M> make_mupairproduction(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, const nlohmann::json& config){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for mupairproduction");
    std::string param_name = config["parametrization"];

    return make_mupairproduction(p_def, medium, cuts, interpol, param_name);
}

} // namespace crosssection
} // namespace PROPOSAL

#undef MUPAIR_PARAM_INTEGRAL_DEC
