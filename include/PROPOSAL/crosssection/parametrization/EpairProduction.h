
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
#include "PROPOSAL/EnergyCutSettings.h"


#define EPAIR_PARAM_INTEGRAL_DEC(param)                                        \
    struct Epair##param : public EpairProductionRhoIntegral {                  \
        Epair##param(bool lpm = false);                                        \
        Epair##param(bool lpm, const ParticleDef&, const Medium&,              \
            double density_correction = 1.0);                                  \
        using base_param_t = EpairProduction;                                  \
                                                                               \
        double FunctionToIntegral(const ParticleDef&, const Component&,        \
            double energy, double v, double r) const;                          \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationName<Epair##param> {                     \
        static constexpr char value[36] = "epair"#param;                       \
    };

namespace PROPOSAL {
namespace crosssection {
    class EpairLPM;
    class EpairProduction : public Parametrization<Component> {

    protected:
        std::shared_ptr<EpairLPM> lpm_;

    public:
        EpairProduction(bool lpm = false);
        EpairProduction(bool lpm, const ParticleDef&, const Medium&,
            double density_correction = 1.0);
        virtual ~EpairProduction() = default;

        //using only_stochastic = std::false_type;
        //using component_wise = std::true_type;

        double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
        KinematicLimits GetKinematicLimits(const ParticleDef&,
            const Component&, double energy) const noexcept override;
    };

    template <> struct ParametrizationName<EpairProduction> {
        static constexpr char value[36] = "epairproduction";
    };

    class EpairProductionRhoIntegral : public EpairProduction {
    public:
        EpairProductionRhoIntegral(bool lpm = false);
        EpairProductionRhoIntegral(bool lpm, const ParticleDef&, const Medium&,
            double density_correction = 1.0);
        virtual ~EpairProductionRhoIntegral() = default;

        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double energy, double v) const override;

        // ----------------------------------------------------------------------------
        /// @brief This is the calculation of the d2Sigma/dvdRo - interface to
        /// Integral
        ///
        /// the function which is defined here is:
        /// \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
        /// \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
        // ----------------------------------------------------------------------------
        virtual double FunctionToIntegral(const ParticleDef&, const Component&,
            double energy, double v, double rho) const = 0;
    };

    EPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)
    EPAIR_PARAM_INTEGRAL_DEC(SandrockSoedingreksoRhode)
    EPAIR_PARAM_INTEGRAL_DEC(ForElectronPositron)

#undef EPAIR_PARAM_INTEGRAL_DEC

    class EpairLPM {
        double mass_;
        double charge_;
        double mol_density_;
        double density_correction_;
        double eLpm_;
        size_t hash;
    public:
        EpairLPM(
            const ParticleDef&, const Medium&, double density_correction = 1.0);
        double suppression_factor(
            double E, double v, double r2, double beta, double xi) const;
        size_t GetHash() const noexcept { return hash;}

    };

    // Factory pattern functions

    /*
    template <typename P, typename M>
    using epair_func_ptr = cross_t_ptr<P, M> (*)(
        P, M, std::shared_ptr<const EnergyCutSettings>, bool, bool, double);

    template <typename Param, typename P, typename M>
    cross_t_ptr<P, M> create_epair(P p_def, M medium,
        std::shared_ptr<const EnergyCutSettings> cuts, bool lpm, bool interpol,
        double density_correction = 1.0)
    {
        auto param = Param(lpm, p_def, medium, density_correction);
        return make_crosssection(param, p_def, medium, cuts, interpol);
    }

    template <typename P, typename M>
    static std::map<std::string, epair_func_ptr<P, M>> epair_map = {
        { "kelnerkokoulinpetrukhin",
            create_epair<EpairKelnerKokoulinPetrukhin, P, M> },
        { "sandrocksoedingreksorhode",
            create_epair<EpairSandrockSoedingreksoRhode, P, M> },
        { "forelectronpositron", create_epair<EpairForElectronPositron, P, M> },
    };

    template <typename P, typename M>
    cross_t_ptr<P, M> make_epairproduction(P p_def, M medium,
        std::shared_ptr<const EnergyCutSettings> cuts, bool interpol, bool lpm,
        const std::string& param_name, double density_correction = 1.0)
    {
        std::string name = param_name;
        std::transform(
            param_name.begin(), param_name.end(), name.begin(), ::tolower);
        auto it = epair_map<P, M>.find(name);
        if (it == epair_map<P, M>.end())
            throw std::logic_error(
                "Unknown parametrization for epairproduction");

        return it->second(
            p_def, medium, cuts, lpm, interpol, density_correction);
    }

    template <typename P, typename M>
    cross_t_ptr<P, M> make_epairproduction(P p_def, M medium,
        std::shared_ptr<const EnergyCutSettings> cuts, bool interpol,
        const nlohmann::json& config, double density_correction = 1.0)
    {
        if (!config.contains("parametrization"))
            throw std::logic_error(
                "No parametrization passed for epairproduction");

        std::string param_name = config["parametrization"];
        bool lpm = config.value("lpm", true);

        return make_epairproduction(
            p_def, medium, cuts, interpol, lpm, param_name, density_correction);
    }
    */

} // namespace crosssection

namespace detail {
    double integrate_dedx_epair(Integral& integral,
        crosssection::EpairProduction const& param, const ParticleDef& p_def,
        const Component& comp, double energy, double v_min, double v_max);

    template <typename T1>
    auto define_epair_dedx_integral(T1 param, ParticleDef const& p_def,
        Component const& comp, EnergyCutSettings const& cut)
    {
        return [param, &p_def, &comp, &cut](Integral& i, double E) {
            auto lim = param.GetKinematicLimits(p_def, comp, E);
            auto v_cut = cut.GetCut(lim, E);
            return integrate_dedx_epair(i, param, p_def, comp, E,
                lim.v_min, v_cut);
        };
    }

    inline auto define_dedx_integral(
        crosssection::EpairSandrockSoedingreksoRhode param,
        ParticleDef const& p_def, Component const& comp,
        EnergyCutSettings const& cut)
    {
        return define_epair_dedx_integral(param, p_def, comp, cut);
    }
    inline auto define_dedx_integral(crosssection::EpairKelnerKokoulinPetrukhin param,
        ParticleDef const& p_def, Component const& comp,
        EnergyCutSettings const& cut)
    {
        return define_epair_dedx_integral(param, p_def, comp, cut);
    }
    inline auto define_dedx_integral(crosssection::EpairForElectronPositron param,
        ParticleDef const& p_def, Component const& comp,
        EnergyCutSettings const& cut)
    {
        return define_epair_dedx_integral(param, p_def, comp, cut);
    }
}
} // namespace PROPOSAL
