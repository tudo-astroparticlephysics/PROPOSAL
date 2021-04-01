
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

#define EPAIR_PARAM_INTEGRAL_DEC(param)                                        \
    struct Epair##param : public EpairProductionRhoIntegral {                  \
        Epair##param(bool lpm = false);                                        \
        Epair##param(bool lpm, const ParticleDef&, const Medium&,              \
            double density_correction = 1.0);                                  \
        using base_param_t = EpairProduction;                                  \
                                                                               \
        std::unique_ptr<Parametrization<Component>> clone() const final;       \
                                                                               \
        double FunctionToIntegral(const ParticleDef&, const Component&,        \
            double energy, double v, double r) const;                          \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationName<Epair##param> {                     \
        static constexpr auto value = #param;                                  \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationId<Epair##param> {                       \
        static constexpr size_t value = 1000000004;                            \
    };

namespace PROPOSAL {
class Integral;
class EnergyCutSettings;
} // namespace PROPOSAL

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


        double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
        KinematicLimits GetKinematicLimits(const ParticleDef&, const Component&,
            double energy) const noexcept override;
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
        size_t GetHash() const noexcept { return hash; }
    };

} // namespace crosssection

} // namespace PROPOSAL
