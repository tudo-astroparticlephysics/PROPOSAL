
#pragma once

#include <cmath>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace PROPOSAL {

/******************************************************************************
*                                   HardBB                                    *
******************************************************************************/

class Interpolant;

class RealPhoton
{
    public:
        RealPhoton() {}
        RealPhoton(const RealPhoton&) {}
        virtual ~RealPhoton() {}

        virtual RealPhoton* clone() const = 0;

        virtual double CalculateHardBB(double energy, double v) = 0;
};

class SoftBB: public RealPhoton
{
    public:
        SoftBB() {}
        SoftBB(const SoftBB&) {}
        virtual ~SoftBB() {}

        RealPhoton* clone() const { return new SoftBB(*this); }

        virtual double CalculateHardBB(double energy, double v);
};

class HardBB: public RealPhoton
{
    public:
        HardBB(const ParticleDef&);
        HardBB(const HardBB&);
        virtual ~HardBB();

        RealPhoton* clone() const { return new HardBB(*this); }

        double CalculateHardBB(double energy, double v);

    private:
        static std::vector<double> x;
        std::vector<Interpolant*> interpolant_;

};

/******************************************************************************
*                                ShadowEffect                                 *
******************************************************************************/

class ShadowEffect
{
    public:
        ShadowEffect() {}
        ShadowEffect(const ShadowEffect&) {}
        virtual ~ShadowEffect() {}

        virtual ShadowEffect* clone() const = 0;

        virtual double CalculateShadowEffect(const Components::Component&, double x, double nu) = 0;
};

class ShadowDutta: public ShadowEffect
{
    public:
        ShadowDutta() {}
        ShadowDutta(const ShadowDutta&) {}
        virtual ~ShadowDutta() {}

        ShadowEffect* clone() const { return new ShadowDutta(*this); }

        double CalculateShadowEffect(const Components::Component&, double x, double nu);
};

class ShadowButkevichMikhailov: public ShadowEffect
{
    public:
        ShadowButkevichMikhailov() {}
        ShadowButkevichMikhailov(const ShadowDutta&) {}
        virtual ~ShadowButkevichMikhailov() {}

        ShadowEffect* clone() const { return new ShadowButkevichMikhailov(*this); }

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

    virtual double DifferentialCrossSection(double energy, double v) = 0;

    virtual double FunctionToDEdxIntegral(double energy, double v);
    virtual double FunctionToDNdxIntegral(double energy, double v);

    virtual IntegralLimits GetIntegralLimits(double energy);
};

} /* PROPOSAL */
