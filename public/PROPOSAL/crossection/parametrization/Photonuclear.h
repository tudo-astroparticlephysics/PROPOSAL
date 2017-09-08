
#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Interpolant.h"

// #define BREMSSTRAHLUNG_DEF(param)                                                                                      \
//     class Brems##param : public Bremsstrahlung                                                                         \
//     {                                                                                                                  \
//         public:                                                                                                        \
//         Brems##param(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());          \
//         Brems##param(const Brems##param&);                                                                             \
//         ~Brems##param();                                                                                               \
//                                                                                                                        \
//         Parametrization* clone() const { return new Brems##param(*this); }                                             \
//                                                                                                                        \
//         double CalculateParametrization(double energy, double v);                                                      \
//     };

namespace PROPOSAL {

/******************************************************************************
*                                   HardBB                                    *
******************************************************************************/

class HardBB
{
    public:
        HardBB(const ParticleDef&);
        HardBB(const HardBB&);
        virtual ~HardBB();

        double CalculateHardBB(double energy, double v);

    private:
        static std::vector<const double> x;

        const ParticleDef particle_def_;
        std::vector<Interpolant*> interpolant_;

};

/******************************************************************************
*                                ShadowEffect                                 *
******************************************************************************/

class ShadowEffect
{
    public:
        ShadowEffect() {}
        virtual ~ShadowEffect() {}

        virtual double CalculateShadowEffect(const Components::Component&, double x, double nu) = 0;
};

class ShadowDutta: public ShadowEffect
{
    public:
        ShadowDutta() {}
        virtual ~ShadowDutta() {}

        double CalculateShadowEffect(const Components::Component&, double x, double nu);
};

class ShadowButkevichMikhailov: public ShadowEffect
{
    public:
        ShadowButkevichMikhailov() {}
        virtual ~ShadowButkevichMikhailov() {}

        double CalculateShadowEffect(const Components::Component&, double x, double nu);
};

/******************************************************************************
*                               Photonuclear                                  *
******************************************************************************/

class Photonuclear : public Parametrization
{
    public:
    Photonuclear(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    Photonuclear(const Photonuclear&);
    virtual ~Photonuclear();

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

    // ----------------------------------------------------------------- //
    // Protected member
    // ----------------------------------------------------------------- //

};

// ------------------------------------------------------------------------- //
// Declare the specific parametrizations
// ------------------------------------------------------------------------- //

// BREMSSTRAHLUNG_DEF(PetrukhinShestakov)
// BREMSSTRAHLUNG_DEF(KelnerKokoulinPetrukhin)
// BREMSSTRAHLUNG_DEF(CompleteScreening)
// BREMSSTRAHLUNG_DEF(AndreevBezrukovBugaev)
// #undef BREMSSTRAHLUNG_DEF

} /* PROPOSAL */
