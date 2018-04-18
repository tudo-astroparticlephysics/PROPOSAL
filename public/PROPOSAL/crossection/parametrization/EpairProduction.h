
#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/methods.h"

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

    // ----------------------------------------------------------------------------
    /// @brief This is the calculation of the d2Sigma/dvdRo - interface to Integral
    ///
    /// the function which is defined here is:
    /// \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    /// \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    // ----------------------------------------------------------------------------
    double FunctionToIntegral(double energy, double r);

    // ----------------------------------------------------------------------------
    /// @brief Landau Pomeranchuk Migdal effect
    ///
    /// Landau Pomeranchuk Migdal effect evaluation,
    /// if the Landau Pomeranchuk Migdal effect is considered in
    /// the calculation, function is modified by a factor
    /// \f[lpm=return=\frac{(1+b)(A+(1+r^2)B)+b(C+(1+r^2)D)+(1-r^2)E}{[(2+r^2)(1+b)
    /// +x(3+r^2)]\ln\Big(1+\frac{1}{x}\Big)+\frac{1-r^2-b}{1+x}-(3+r^2)}\f]
    // ----------------------------------------------------------------------------
    double lpm(double energy, double r2, double b, double x);

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    const std::string& GetName() const { return name_; }
    virtual size_t GetHash() const;

protected:
    bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const;

    static const std::string name_;

    double v_;
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

    Parametrization* clone() const { return new EpairProductionRhoIntegral(*this); }

    virtual double DifferentialCrossSection(double energy, double v);

private:
    bool compare(const Parametrization&) const;

    Integral integral_;
};

class EpairProductionRhoInterpolant : public EpairProductionRhoIntegral
{
public:
    EpairProductionRhoInterpolant(const ParticleDef&,
                                  const Medium&,
                                  const EnergyCutSettings&,
                                  double multiplier,
                                  bool lpm,
                                  InterpolationDef = InterpolationDef());
    EpairProductionRhoInterpolant(const EpairProductionRhoInterpolant&);
    virtual ~EpairProductionRhoInterpolant();

    Parametrization* clone() const { return new EpairProductionRhoInterpolant(*this); }

    double DifferentialCrossSection(double energy, double v);

private:
    bool compare(const Parametrization&) const;

    double FunctionToBuildEpairInterpolant(double energy, double v, int component);

    std::vector<Interpolant*> interpolant_; //!< Interpolates function used by dNdx and dEdx
};

} // namespace PROPOSAL
