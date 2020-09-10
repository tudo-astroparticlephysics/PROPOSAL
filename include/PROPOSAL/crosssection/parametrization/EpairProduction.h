
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
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/crosssection/CrossSection.h"

using PROPOSAL::Components::Component;

#define EPAIR_PARAM_INTEGRAL_DEC(param)                                        \
    struct Epair##param : public EpairProductionRhoIntegral {                  \
        Epair##param(bool lpm);                                                \
        using base_param_t = EpairProduction;                                  \
                                                                               \
        double FunctionToIntegral(const ParticleDef&, const Component&,        \
            double energy, double v, double lpm) const;                              \
    };

namespace PROPOSAL {
namespace crosssection {
class EpairProduction : public Parametrization {

protected:
    // ----------------------------------------------------------------------------
    /// @brief Landau Pomeranchuk Migdal effect
    ///
    /// Landau Pomeranchuk Migdal effect evaluation,
    /// if the Landau Pomeranchuk Migdal effect is considered in
    /// the calculation, function is modified by a factor
    /// \f[lpm=return=\frac{(1+b)(A+(1+r^2)B)+b(C+(1+r^2)D)+(1-r^2)E}{[(2+r^2)(1+b)
    /// +x(3+r^2)]\ln\Big(1+\frac{1}{x}\Big)+\frac{1-r^2-b}{1+x}-(3+r^2)}\f]
    // ----------------------------------------------------------------------------
    double lpm(const ParticleDef&, const Medium&, double energy, double v,
        double r2, double b, double x);

    bool init_lpm_effect_;
    bool lpm_;
    double eLpm_;

public:
    EpairProduction(bool);
    virtual ~EpairProduction() = default;

    using only_stochastic = std::false_type;
    using component_wise = std::true_type;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double energy) const noexcept override;
    size_t GetHash() const noexcept override;
};

class EpairProductionRhoIntegral : public EpairProduction {
    public:
    EpairProductionRhoIntegral(bool);
    virtual ~EpairProductionRhoIntegral() = default;

    double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v) const override;

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the d2Sigma/dvdRo - interface to
    /// Integral
    ///
    /// the function which is defined here is:
    /// \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    /// \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    // ----------------------------------------------------------------------------
    virtual double FunctionToIntegral(const ParticleDef&, const Component&,
        double energy, double v, double rho) const
        = 0;

};

EPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)
EPAIR_PARAM_INTEGRAL_DEC(SandrockSoedingreksoRhode)

#undef EPAIR_PARAM_INTEGRAL_DEC

// Factory pattern functions

template <typename P, typename M>
using epair_func_ptr = cross_t_ptr<P, M>(*)(P, M, std::shared_ptr<const
        EnergyCutSettings>, bool, bool);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_epair(P p_def, M medium,std::shared_ptr<const
        EnergyCutSettings> cuts, bool lpm, bool interpol) {
    auto param = Param(lpm);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template<typename P, typename M>
static std::map<std::string, epair_func_ptr<P, M>> epair_map = {
        {"KelnerKokoulinPetrukhin", create_epair<EpairKelnerKokoulinPetrukhin, P, M>},
        {"SandrockSoedingreksoRhode", create_epair<EpairSandrockSoedingreksoRhode, P, M>},
};

template<typename P, typename M>
cross_t_ptr<P, M> make_epairproduction(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, const nlohmann::json& config){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for epairproduction");

    std::string param_name = config["parametrization"];
    auto it = epair_map<P, M>.find(param_name);
    if (it == epair_map<P, M>.end())
        throw std::logic_error("Unknown parametrization for epairproduction");

    bool lpm = config.value("lpm", true);
    return it->second(p_def, medium, cuts, lpm, interpol);
}

} // namespace crosssection

template <>
double integrate_dedx(Integral&, crosssection::EpairProduction&, const ParticleDef&,
    const Component&, double, double, double);

} // namespace PROPOSAL
