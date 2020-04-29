
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

    virtual double CalculateHardComponent(double energy, double v) = 0;

    virtual const std::string& GetName() const = 0;

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
    virtual ~ShadowEffect() {}


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
    Photonuclear(const ParticleDef&, const component_list&);

    virtual double DifferentialCrossSection(double energy, double v) = 0;

    virtual KinematicLimits GetKinematicLimits(double energy);

protected:
};

} // namespace PROPOSAL
