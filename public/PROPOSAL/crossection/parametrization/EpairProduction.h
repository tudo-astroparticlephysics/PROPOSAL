
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

#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"

#define EPAIR_PARAM_INTEGRAL_DEC(param)                                                                                \
    class Epair##param : public EpairProductionRhoIntegral                                                             \
    {                                                                                                                  \
    public:                                                                                                            \
        Epair##param(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier, bool lpm);                  \
        Epair##param(const Epair##param&);                                                                             \
        virtual ~Epair##param();                                                                                       \
                                                                                                                       \
        virtual Parametrization* clone() const { return new Epair##param(*this); }                                     \
        static EpairProduction* create(const ParticleDef& particle_def,                                                \
                                       std::shared_ptr<const Medium> medium,                                           \
                                       double multiplier,                                                              \
                                       bool lpm)                                                                       \
        {                                                                                                              \
            return new Epair##param(particle_def, medium, multiplier, lpm);                                            \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToIntegral(double energy, double v, double lpm);                                                \
                                                                                                                       \
        const std::string& GetName() const { return name_; }                                                           \
                                                                                                                       \
    protected:                                                                                                         \
        static const std::string name_;                                                                                \
    };


namespace PROPOSAL {

class Interpolant;

class EpairProduction : public Parametrization
{
public:
    EpairProduction(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier, bool lpm);
    EpairProduction(const EpairProduction&);
    virtual ~EpairProduction();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the dSigma/dv
    ///
    /// \f[e_{Pair}=return = \rho N_Z z^2 \Big[ \int_{1-r_{max}}^{aux}
    /// f(r)dr + \int^1_{aux} f(r) dr \Big]\f]
    /// with \f$ aux=max(1-r_{Max}, ComputerPrecision)\f$ and \f$r_{max} =
    /// \sqrt{1-\frac{4m_e}{E_p v}}\Big(1-\frac{6m_p^2}{E_p^2(1-v)}\Big)\f$
    ///
    // ----------------------------------------------------------------------------
    virtual InteractionType GetInteractionType() const final {return InteractionType::Epair;}
    virtual double DifferentialCrossSection(double energy, double v) = 0;

    virtual KinematicLimits GetKinematicLimits(double energy);


protected:
    bool compare(const Parametrization&) const;

    // ----------------------------------------------------------------------------
    /// @brief Landau Pomeranchuk Migdal effect
    ///
    /// Landau Pomeranchuk Migdal effect evaluation,
    /// if the Landau Pomeranchuk Migdal effect is considered in
    /// the calculation, function is modified by a factor
    /// \f[lpm=return=\frac{(1+b)(A+(1+r^2)B)+b(C+(1+r^2)D)+(1-r^2)E}{[(2+r^2)(1+b)
    /// +x(3+r^2)]\ln\Big(1+\frac{1}{x}\Big)+\frac{1-r^2-b}{1+x}-(3+r^2)}\f]
    // ----------------------------------------------------------------------------
    double lpm(double energy, double v, double r2, double b, double x);

    bool init_lpm_effect_;
    bool lpm_;
    double eLpm_;
};

// ------------------------------------------------------------------------- //
// Differentiate between rho integration & interpolation
// ------------------------------------------------------------------------- //

class EpairProductionRhoIntegral : public EpairProduction
{
public:
    EpairProductionRhoIntegral(const ParticleDef&,
                               std::shared_ptr<const Medium>,
                               double multiplier,
                               bool lpm);
    EpairProductionRhoIntegral(const EpairProductionRhoIntegral&);
    virtual ~EpairProductionRhoIntegral();

    Parametrization* clone() const = 0;

    virtual double DifferentialCrossSection(double energy, double v);

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the d2Sigma/dvdRo - interface to Integral
    ///
    /// the function which is defined here is:
    /// \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    /// \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    // ----------------------------------------------------------------------------
    virtual double FunctionToIntegral(double energy, double v, double rho) = 0;

    virtual size_t GetHash() const;

private:
    bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const;

    Integral integral_;
};

/******************************************************************************
 *                     Declare Integral Parametrizations                      *
 ******************************************************************************/

EPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)
EPAIR_PARAM_INTEGRAL_DEC(SandrockSoedingreksoRhode)

#undef EPAIR_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
