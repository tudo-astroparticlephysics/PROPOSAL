
#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace PROPOSAL {
class Ionization : public Parametrization
{
    public:
    Ionization(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    Ionization(const Ionization&);
    virtual ~Ionization();

    virtual Parametrization* clone() const { return new Ionization(*this); }

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    double DifferentialCrossSection(double energy, double v);
    double CalculateParametrization(double energy, double v);

    double FunctionToDEdxIntegral(double energy, double v);
    double FunctionToDNdxIntegral(double energy, double v);

    IntegralLimits GetIntegralLimits(double energy);

    double InelCorrection(double energy, double v);
};


} /* PROPOSAL */
