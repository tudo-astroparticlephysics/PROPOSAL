
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

        virtual const std::string& GetName() const = 0;
};

class SoftBB: public RealPhoton
{
    public:
        SoftBB();
        SoftBB(const SoftBB& bb);
        virtual ~SoftBB();

        RealPhoton* clone() const { return new SoftBB(*this); }

        virtual double CalculateHardBB(double energy, double v);

        virtual const std::string& GetName() const { return name_; }

    private:
        static const std::string name_;
};

class HardBB: public RealPhoton
{
    public:
        HardBB(const ParticleDef&);
        HardBB(const HardBB&);
        virtual ~HardBB();

        RealPhoton* clone() const { return new HardBB(*this); }

        double CalculateHardBB(double energy, double v);

        virtual const std::string& GetName() const { return name_; }

    private:
        static std::vector<double> x;
        std::vector<Interpolant*> interpolant_;

        static const std::string name_;

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

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    virtual const std::string& GetName() const = 0;
    virtual size_t GetHash() const = 0;
};

class ShadowDuttaRenoSarcevicSeckel: public ShadowEffect
{
    public:
    ShadowDuttaRenoSarcevicSeckel() {}
    ShadowDuttaRenoSarcevicSeckel(const ShadowDuttaRenoSarcevicSeckel&) {}
    virtual ~ShadowDuttaRenoSarcevicSeckel() {}

    ShadowEffect* clone() const { return new ShadowDuttaRenoSarcevicSeckel(*this); }
    static ShadowEffect* create() { return new ShadowDuttaRenoSarcevicSeckel(); }

    double CalculateShadowEffect(const Components::Component&, double x, double nu);

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    virtual const std::string& GetName() const { return name_; }
    virtual size_t GetHash() const;

    private:
    static const std::string name_;
};

class ShadowButkevichMikhailov: public ShadowEffect
{
    public:
    ShadowButkevichMikhailov() {}
    ShadowButkevichMikhailov(const ShadowButkevichMikhailov&) {}
    virtual ~ShadowButkevichMikhailov() {}

    ShadowEffect* clone() const { return new ShadowButkevichMikhailov(*this); }
    static ShadowEffect* create() { return new ShadowButkevichMikhailov(); }

    double CalculateShadowEffect(const Components::Component&, double x, double nu);

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    virtual const std::string& GetName() const { return name_; }
    virtual size_t GetHash() const;

    private:
    static const std::string name_;
};

/******************************************************************************
*                               Photonuclear                                  *
******************************************************************************/

class Photonuclear : public Parametrization
{
    public:
    Photonuclear(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier);
    Photonuclear(const Photonuclear&);
    virtual ~Photonuclear();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v) = 0;

    virtual IntegralLimits GetIntegralLimits(double energy);

    protected:
    virtual bool compare(const Parametrization&) const;

};

} /* PROPOSAL */
