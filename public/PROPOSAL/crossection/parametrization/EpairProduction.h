
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

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"

using PROPOSAL::Components::Component;

#define EPAIR_PARAM_INTEGRAL_DEC(param)                                        \
    struct Epair##param : public EpairProductionRhoIntegral {                  \
        Epair##param(bool lpm);                                                \
        using base_param_t = EpairProduction;                                  \
                                                                               \
        double FunctionToIntegral(const ParticleDef&, const Component&,        \
            double energy, double v, double lpm);                              \
    };

namespace PROPOSAL {

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

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the dSigma/dv
    ///
    /// \f[e_{Pair}=return = \rho N_Z z^2 \Big[ \int_{1-r_{max}}^{aux}
    /// f(r)dr + \int^1_{aux} f(r) dr \Big]\f]
    /// with \f$ aux=max(1-r_{Max}, ComputerPrecision)\f$ and \f$r_{max} =
    /// \sqrt{1-\frac{4m_e}{E_p v}}\Big(1-\frac{6m_p^2}{E_p^2(1-v)}\Big)\f$
    ///
    // ----------------------------------------------------------------------------
    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v)
        = 0;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double energy) const noexcept override;
    size_t GetHash() const noexcept override;
};

template <>
double integrate_dedx(Integral&, EpairProduction&, const ParticleDef&,
    const Component&, double, double, double);

class EpairProductionRhoIntegral : public EpairProduction {
private:
    Integral integral_;

public:
    EpairProductionRhoIntegral(bool);
    virtual ~EpairProductionRhoIntegral() = default;

    double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v) override;

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the d2Sigma/dvdRo - interface to
    /// Integral
    ///
    /// the function which is defined here is:
    /// \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    /// \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    // ----------------------------------------------------------------------------
    virtual double FunctionToIntegral(const ParticleDef&, const Component&,
        double energy, double v, double rho)
        = 0;

};

EPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)
EPAIR_PARAM_INTEGRAL_DEC(SandrockSoedingreksoRhode)

#undef EPAIR_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
