
#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#define BREMSSTRAHLUNG_DEF(param)                                                                                      \
    class Brems##param : public Bremsstrahlung                                                                         \
    {                                                                                                                  \
        public:                                                                                                        \
        Brems##param(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());          \
        Brems##param(const Brems##param&);                                                                             \
        ~Brems##param();                                                                                               \
                                                                                                                       \
        Parametrization* clone() const { return new Brems##param(*this); }                                             \
        static Parametrization* create(const ParticleDef& particle_def,                                                \
                                       const Medium& medium,                                                           \
                                       const EnergyCutSettings& cuts,                                                  \
                                       Definition def = Definition())                                                  \
        {                                                                                                              \
            return new Brems##param(particle_def, medium, cuts, def);                                                  \
        }                                                                                                              \
                                                                                                                       \
        double CalculateParametrization(double energy, double v);                                                      \
    };

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

    double lpm(double energy, double v);

    // ----------------------------------------------------------------- //
    // Protected member
    // ----------------------------------------------------------------- //

    bool lorenz_;       /// enable lorenz cut
    double lorenz_cut_; /// in [MeV] // - set to 1.e6 in Constructor
    double eLpm_;
};

// ------------------------------------------------------------------------- //
// Declare the specific parametrizations
// ------------------------------------------------------------------------- //

BREMSSTRAHLUNG_DEF(PetrukhinShestakov)
BREMSSTRAHLUNG_DEF(KelnerKokoulinPetrukhin)
BREMSSTRAHLUNG_DEF(CompleteScreening)
BREMSSTRAHLUNG_DEF(AndreevBezrukovBugaev)

#undef BREMSSTRAHLUNG_DEF

} /* PROPOSAL */
