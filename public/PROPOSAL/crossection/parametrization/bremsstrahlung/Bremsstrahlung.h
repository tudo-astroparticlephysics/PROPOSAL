
#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace PROPOSAL {
class Bremsstrahlung : public Parametrization
{
    public:
    Bremsstrahlung(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    Bremsstrahlung(const Bremsstrahlung&);
    virtual ~Bremsstrahlung();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v);
    virtual double CalculateParametrization(double energy, double v) = 0;

    virtual double FunctionToDEdxIntegral(double energy, double v);
    virtual double FunctionToDNdxIntegral(double energy, double v);
    virtual IntegralLimits GetIntegralLimits(double energy);

    protected:

    // ----------------------------------------------------------------- //
    // Protected methods
    // ----------------------------------------------------------------- //

    double lpm(double energy, double v, double s1);

    // ----------------------------------------------------------------- //
    // Protected member
    // ----------------------------------------------------------------- //

    bool lorenz_;       /// enable lorenz cut
    double lorenz_cut_; /// in [MeV] // - set to 1.e6 in Constructor
    double eLpm_;
};
} /* PROPOSAL */
