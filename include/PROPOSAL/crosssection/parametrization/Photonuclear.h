
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
#include <memory>

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

using std::shared_ptr;

namespace PROPOSAL {
class Interpolant;
namespace crosssection {
class RealPhoton {
public:
    RealPhoton() = default;
    virtual ~RealPhoton() {}

    virtual double CalculateHardComponent(double energy, double v) = 0;

    virtual const std::string& GetName() const = 0;
};

class SoftComponent : public RealPhoton {
    static const std::string name_;

public:
    SoftComponent() = default;
    virtual ~SoftComponent() = default;

    virtual double CalculateHardComponent(double energy, double v);

    virtual const std::string& GetName() const { return name_; }
};

class HardComponent : public RealPhoton {
    static std::vector<double> x;
    std::vector<shared_ptr<Interpolant>> interpolant_;

    static const std::string name_;

public:
    HardComponent(const ParticleDef&);
    virtual ~HardComponent() = default;

    double CalculateHardComponent(double energy, double v);

    virtual const std::string& GetName() const { return name_; }
};

class ShadowEffect {
public:
    ShadowEffect() {}
    virtual ~ShadowEffect() {}

    virtual double CalculateShadowEffect(const Component&, double x, double nu)
        = 0;

    virtual const std::string& GetName() const = 0;
    virtual size_t GetHash() const = 0;
};

class ShadowDuttaRenoSarcevicSeckel : public ShadowEffect {
public:
    ShadowDuttaRenoSarcevicSeckel()
        : ShadowEffect()
    {
    }

    double CalculateShadowEffect(const Component&, double x, double nu);

    virtual const std::string& GetName() const { return name_; }
    virtual size_t GetHash() const;

private:
    static const std::string name_;
};

class ShadowButkevichMikhailov : public ShadowEffect {
public:
    ShadowButkevichMikhailov()
        : ShadowEffect()
    {
    }

    double CalculateShadowEffect(const Component&, double x, double nu);

    virtual const std::string& GetName() const { return name_; }
    virtual size_t GetHash() const;

private:
    static const std::string name_;
};

class Photonuclear : public Parametrization {
public:
    Photonuclear();
    virtual ~Photonuclear() = default;

    using only_stochastic = std::false_type;
    using component_wise = std::true_type;

    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double, double) const
        = 0;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double) const noexcept override;
};

} // namespace crosssection
} // namespace PROPOSAL
