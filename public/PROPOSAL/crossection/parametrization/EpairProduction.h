
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

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

#define EPAIR_PARAM_INTEGRAL_DEC(param)                                                                                \
    class Epair##param : public EpairProductionRhoIntegral                                                             \
    {                                                                                                                  \
    public:                                                                                                            \
        Epair##param(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, bool lpm);        \
        Epair##param(const Epair##param&);                                                                             \
        virtual ~Epair##param();                                                                                       \
                                                                                                                       \
        virtual Parametrization* clone() const { return new Epair##param(*this); }                                     \
        static EpairProduction* create(const ParticleDef& particle_def,                                                \
                                       const Medium& medium,                                                           \
                                       const EnergyCutSettings& cuts,                                                  \
                                       double multiplier,                                                              \
                                       bool lpm)                                                                       \
        {                                                                                                              \
            return new Epair##param(particle_def, medium, cuts, multiplier, lpm);                                      \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToIntegral(double energy, double v, double lpm);                                              \
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
    EpairProduction(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, bool lpm);
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
    virtual double DifferentialCrossSection(double energy, double v) = 0;

    virtual IntegralLimits GetIntegralLimits(double energy);


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
                               const Medium&,
                               const EnergyCutSettings&,
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

/******************************************************************************
 *                    Declare Interpolant Parametrizations                    *
 ******************************************************************************/

template<class Param = EpairKelnerKokoulinPetrukhin>
class EpairProductionRhoInterpolant : public Param
{
public:
    typedef std::vector<Interpolant*> InterpolantVec;

public:
    EpairProductionRhoInterpolant(const ParticleDef&,
                       const Medium&,
                       const EnergyCutSettings&,
                       double multiplier,
                       bool lpm,
                       InterpolationDef def = InterpolationDef());
    EpairProductionRhoInterpolant(const EpairProductionRhoInterpolant&);
    virtual ~EpairProductionRhoInterpolant();

    Parametrization* clone() const { return new EpairProductionRhoInterpolant<Param>(*this); }
    static EpairProduction* create(const ParticleDef& particle_def,
                                const Medium& medium,
                                const EnergyCutSettings& cuts,
                                double multiplier,
                                bool lpm,
                                InterpolationDef def = InterpolationDef())
    {
        return new EpairProductionRhoInterpolant<Param>(particle_def, medium, cuts, multiplier, lpm, def);
    }

    double DifferentialCrossSection(double energy, double v);

protected:
    virtual bool compare(const Parametrization&) const;
    double FunctionToBuildPhotoInterpolant(double energy, double v, int component);

    InterpolantVec interpolant_;
};

template<class Param>
EpairProductionRhoInterpolant<Param>::EpairProductionRhoInterpolant(const ParticleDef& particle_def,
                                              const Medium& medium,
                                              const EnergyCutSettings& cuts,
                                              double multiplier,
                                              bool lpm,
                                              InterpolationDef def)
    : Param(particle_def, medium, cuts, multiplier, lpm)
    , interpolant_(this->medium_->GetNumComponents(), NULL)
{
    std::vector<Interpolant2DBuilder> builder2d(this->components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(this->components_.size());

    for (unsigned int i = 0; i < this->components_.size(); ++i)
    {
        builder2d[i]
            .SetMax1(def.nodes_cross_section)
            .SetX1Min(this->particle_def_.low)
            .SetX1Max(def.max_node_energy)
            .SetMax2(def.nodes_cross_section)
            .SetX2Min(0.0)
            .SetX2Max(1.0)
            .SetRomberg1(def.order_of_interpolation)
            .SetRational1(false)
            .SetRelative1(false)
            .SetIsLog1(true)
            .SetRomberg2(def.order_of_interpolation)
            .SetRational2(false)
            .SetRelative2(false)
            .SetIsLog2(false)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction2D(std::bind(&EpairProductionRhoInterpolant::FunctionToBuildPhotoInterpolant, this, std::placeholders::_1, std::placeholders::_2, i));

        builder_container2d[i].first  = &builder2d[i];
        builder_container2d[i].second = &interpolant_[i];
    }

    Helper::InitializeInterpolation("Epair", builder_container2d, std::vector<Parametrization*>(1, this), def);
}

template<class Param>
EpairProductionRhoInterpolant<Param>::EpairProductionRhoInterpolant(const EpairProductionRhoInterpolant& photo)
    : Param(photo)
    , interpolant_()
{
    interpolant_.resize(photo.interpolant_.size());

    for (unsigned int i = 0; i < photo.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*photo.interpolant_[i]);
    }
}

template<class Param>
EpairProductionRhoInterpolant<Param>::~EpairProductionRhoInterpolant()
{
    for (std::vector<Interpolant*>::const_iterator iter = interpolant_.begin(); iter != interpolant_.end(); ++iter)
    {
        delete *iter;
    }

    interpolant_.clear();
}

template<class Param>
bool EpairProductionRhoInterpolant<Param>::compare(const Parametrization& parametrization) const
{
    const EpairProductionRhoInterpolant<Param>* epair = static_cast<const EpairProductionRhoInterpolant<Param>*>(&parametrization);

    if (interpolant_.size() != epair->interpolant_.size())
        return false;

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *epair->interpolant_[i])
            return false;
    }

    return EpairProduction::compare(parametrization);
}

template<class Param>
double EpairProductionRhoInterpolant<Param>::DifferentialCrossSection(double energy, double v)
{
    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);

    if (v >= limits.vUp)
    {
        return std::max(
            interpolant_.at(this->component_index_)->Interpolate(energy, std::log(v / limits.vUp) / std::log(limits.vMax / limits.vUp)),
            0.0);
    } else
    {
        return Param::DifferentialCrossSection(energy, v);
    }
}

template<class Param>
double EpairProductionRhoInterpolant<Param>::FunctionToBuildPhotoInterpolant(double energy, double v, int component)
{
    this->component_index_                 = component;
    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * std::exp(v * std::log(limits.vMax / limits.vUp));

    return Param::DifferentialCrossSection(energy, v);
}

#undef EPAIR_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
