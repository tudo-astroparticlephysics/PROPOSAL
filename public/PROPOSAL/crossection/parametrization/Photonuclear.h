
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

#include <cmath>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace PROPOSAL {

/******************************************************************************
 *                               HardComponent                                 *
 ******************************************************************************/

class Interpolant;

class RealPhoton
{
public:
    RealPhoton() {}
    RealPhoton(const RealPhoton&) {}
    virtual ~RealPhoton() {}

    virtual RealPhoton* clone() const = 0;

    bool operator==(const RealPhoton&) const;
    bool operator!=(const RealPhoton&) const;

    virtual double CalculateHardComponent(double energy, double v) = 0;

    virtual const std::string& GetName() const = 0;

protected:
    virtual bool compare(const RealPhoton&) const;
};

class SoftComponent : public RealPhoton
{
public:
    SoftComponent();
    SoftComponent(const SoftComponent&);
    virtual ~SoftComponent();

    RealPhoton* clone() const { return new SoftComponent(*this); }

    virtual double CalculateHardComponent(double energy, double v);

    virtual const std::string& GetName() const { return name_; }

private:
    static const std::string name_;
};

class HardComponent : public RealPhoton
{
public:
    HardComponent(const ParticleDef&);
    HardComponent(const HardComponent&);
    virtual ~HardComponent();

    RealPhoton* clone() const { return new HardComponent(*this); }

    double CalculateHardComponent(double energy, double v);

    virtual const std::string& GetName() const { return name_; }

private:
    virtual bool compare(const RealPhoton&) const;

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

    bool operator==(const ShadowEffect&) const;
    bool operator!=(const ShadowEffect&) const;

    virtual double CalculateShadowEffect(const Components::Component&, double x, double nu) = 0;

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    virtual const std::string& GetName() const = 0;
    virtual size_t GetHash() const             = 0;
};

class ShadowDuttaRenoSarcevicSeckel : public ShadowEffect
{
public:
    ShadowDuttaRenoSarcevicSeckel()
        : ShadowEffect()
    {
    }
    ShadowDuttaRenoSarcevicSeckel(const ShadowDuttaRenoSarcevicSeckel& sh)
        : ShadowEffect(sh)
    {
    }
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

class ShadowButkevichMikhailov : public ShadowEffect
{
public:
    ShadowButkevichMikhailov()
        : ShadowEffect()
    {
    }
    ShadowButkevichMikhailov(const ShadowButkevichMikhailov& sh)
        : ShadowEffect(sh)
    {
    }
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

} // namespace PROPOSAL
