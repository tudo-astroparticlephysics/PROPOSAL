
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

#define MUPAIR_PARAM_INTEGRAL_DEC(param)                                                                               \
    class Mupair##param : public MupairProductionRhoIntegral                                                           \
    {                                                                                                                  \
    public:                                                                                                            \
        Mupair##param(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier);                           \
        Mupair##param(const Mupair##param&);                                                                           \
        virtual ~Mupair##param();                                                                                      \
                                                                                                                       \
        virtual Parametrization* clone() const { return new Mupair##param(*this); }                                    \
        static MupairProduction* create(const ParticleDef& particle_def,                                               \
                                       std::shared_ptr<const Medium> medium,                                           \
                                       double multiplier)                                                              \
        {                                                                                                              \
            return new Mupair##param(particle_def, medium, multiplier);                                                \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToIntegral(double energy, double v, double r);                                                  \
                                                                                                                       \
        const std::string& GetName() const { return name_; }                                                           \
                                                                                                                       \
    protected:                                                                                                         \
        static const std::string name_;                                                                                \
    };


namespace PROPOSAL {

class Interpolant;

class MupairProduction : public Parametrization
{
public:
    MupairProduction(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier);
    MupairProduction(const MupairProduction&);
    virtual ~MupairProduction();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the dSigma/dv
    // ----------------------------------------------------------------------------
    virtual InteractionType GetInteractionType() const final {return InteractionType::MuPair;}
    virtual double DifferentialCrossSection(double energy, double v) = 0;
    virtual double FunctionToIntegral(double energy, double v, double rho) = 0;
    virtual double Calculaterho(double energy, double v, double rnd1, double rnd2);

    virtual KinematicLimits GetKinematicLimits(double energy);

protected:
    bool compare(const Parametrization&) const;
    Integral drho_integral_;
};

// ------------------------------------------------------------------------- //
// Differentiate between rho integration & interpolation
// ------------------------------------------------------------------------- //

class MupairProductionRhoIntegral : public MupairProduction
{
public:
    MupairProductionRhoIntegral(const ParticleDef&,
                               std::shared_ptr<const Medium>,
                               double multiplier);
    MupairProductionRhoIntegral(const MupairProductionRhoIntegral&);
    virtual ~MupairProductionRhoIntegral();

    Parametrization* clone() const = 0;

    virtual double DifferentialCrossSection(double energy, double v);

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the d2Sigma/dvdRo - interface to Integral
    ///
    // ----------------------------------------------------------------------------

    virtual size_t GetHash() const;

private:
    bool compare(const Parametrization&) const;
    //virtual void print(std::ostream&) const;

    Integral integral_;
};

/******************************************************************************
 *                     Declare Integral Parametrizations                      *
 ******************************************************************************/

MUPAIR_PARAM_INTEGRAL_DEC(KelnerKokoulinPetrukhin)

#undef MUPAIR_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
