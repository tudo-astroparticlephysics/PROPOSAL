
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

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"


namespace PROPOSAL {

class PhotoPairProduction : public Parametrization {
public:
    PhotoPairProduction(const ParticleDef&, const component_list&);
    virtual double DifferentialCrossSection(double energy, double x) = 0;

    virtual KinematicLimits GetKinematicLimits(double energy);
};

class PhotoPairTsai : public PhotoPairProduction {
public:
    PhotoPairTsai(const ParticleDef&, const component_list&);

    virtual double DifferentialCrossSection(double energy, double x);
};

class PhotoAngleDistribution {
protected:
    component_list components_;
    Components::Component current_component_;

public:
    PhotoAngleDistribution(const component_list&);

    struct DeflectionAngles {
        double cosphi0, theta0, cosphi1, theta1;
    };

    virtual DeflectionAngles SampleAngles(double energy, double rho) = 0;
};

class PhotoAngleTsaiIntegral : public PhotoAngleDistribution {
    Integral integral_;
    double FunctionToIntegral(double energy, double x, double theta);

public:
    PhotoAngleTsaiIntegral(const component_list&);

    DeflectionAngles SampleAngles(double, double) override;
};

class PhotoAngleNoDeflection : public PhotoAngleDistribution {
public:
    PhotoAngleNoDeflection(const component_list&);

    DeflectionAngles SampleAngles(double, double) override;
};

class PhotoAngleEGS : public PhotoAngleDistribution {
public:
    PhotoAngleEGS(const component_list&);

    DeflectionAngles SampleAngles(double, double) override;
};

} // namespace PROPOSAL
