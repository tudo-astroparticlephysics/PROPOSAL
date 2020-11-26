
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
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"

using PROPOSAL::Components::Component;

#define EPAIR_PARAM_INTEGRAL_DEC(param)                                        \
    struct Epair##param : public EpairProductionRhoIntegral {                  \
        Epair##param(bool lpm = false);                                        \
        Epair##param(bool lpm, const ParticleDef&, const Medium&, double       \
        density_correction = 1.0);                                             \
        using base_param_t = EpairProduction;                                  \
                                                                               \
        double FunctionToIntegral(const ParticleDef&, const Component&,        \
            double energy, double v, double r) const;                          \
    };

namespace PROPOSAL {
namespace crosssection {
class EpairLPM;
class EpairProduction : public Parametrization {

protected:
    std::shared_ptr<EpairLPM> lpm_;

public:
    EpairProduction(bool lpm = false);
    EpairProduction(bool lpm, const ParticleDef&, const Medium&,
                    double density_correction = 1.0);
    virtual ~EpairProduction() = default;

    using only_stochastic = std::false_type;
    using component_wise = std::true_type;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    std::tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double energy) const noexcept override;
    size_t GetHash() const noexcept override;
};

class EpairProductionRhoIntegral : public EpairProduction {
    public:
    EpairProductionRhoIntegral(bool lpm = false);
    EpairProductionRhoIntegral(bool lpm, const ParticleDef&, const Medium&,
                               double density_correction = 1.0);
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
EPAIR_PARAM_INTEGRAL_DEC(ForElectronPositron)

#undef EPAIR_PARAM_INTEGRAL_DEC

class EpairLPM {
public:
    EpairLPM(const ParticleDef&, const Medium&, double density_correction=1.0);
    double suppression_factor(double E, double v, double r2, double beta,
                              double xi) const;
    size_t GetHash() const noexcept;
private:
    double mass_;
    double charge_;
    double mol_density_;
    double density_correction_;
    double eLpm_;
};

// Factory pattern functions

template <typename P, typename M>
using epair_func_ptr = cross_t_ptr<P, M>(*)(P, M, std::shared_ptr<const
        EnergyCutSettings>, bool, bool, double);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_epair(P p_def, M medium,std::shared_ptr<const
        EnergyCutSettings> cuts, bool lpm, bool interpol,
        double density_correction = 1.0) {
    auto param = Param(lpm, p_def, medium, density_correction);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template<typename P, typename M>
static std::map<std::string, epair_func_ptr<P, M>> epair_map = {
        {"kelnerkokoulinpetrukhin", create_epair<EpairKelnerKokoulinPetrukhin, P, M>},
        {"sandrocksoedingreksorhode", create_epair<EpairSandrockSoedingreksoRhode, P, M>},
        {"forelectronpositron", create_epair<EpairForElectronPositron, P, M>},
};

template<typename P, typename M>
cross_t_ptr<P, M> make_epairproduction(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, bool lpm,
        const std::string& param_name, double density_correction = 1.0){
    std::string name = param_name;
    std::transform(param_name.begin(), param_name.end(), name.begin(), ::tolower);
    auto it = epair_map<P, M>.find(name);
    if (it == epair_map<P, M>.end())
        throw std::logic_error("Unknown parametrization for epairproduction");

    return it->second(p_def, medium, cuts, lpm, interpol, density_correction);
}

template<typename P, typename M>
cross_t_ptr<P, M> make_epairproduction(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, const nlohmann::json& config,
        double density_correction = 1.0){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for epairproduction");

    std::string param_name = config["parametrization"];
    bool lpm = config.value("lpm", true);

    return make_epairproduction(p_def, medium, cuts, interpol, lpm, param_name,
                                density_correction);
}

} // namespace crosssection

template <>
double integrate_dedx(Integral&, crosssection::EpairProduction&, const ParticleDef&,
    const Component&, double, double, double);

} // namespace PROPOSAL
