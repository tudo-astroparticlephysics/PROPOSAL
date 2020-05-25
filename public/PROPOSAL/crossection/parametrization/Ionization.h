
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

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {
class Ionization : public Parametrization {
    using only_stochastic = std::false_type;
protected:
    EnergyCutSettings cuts_;

    double Delta(const Medium&, double beta, double gamma);

public:
    Ionization(const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(const ParticleDef&, const Component&, double);
    virtual KinematicLimits GetKinematicLimits( const ParticleDef&, const Medium&, double energy) = 0;
    double DifferentialCrossSection(const ParticleDef&, const Component&, double, double);
};

class IonizBetheBlochRossi : public Ionization {
    double InelCorrection(const ParticleDef&, const Medium&, double, double);
    double CrossSectionWithoutInelasticCorrection(
        const ParticleDef&, const Medium&, double, double);

public:
    IonizBetheBlochRossi(const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(const ParticleDef&, const Medium&, double);
    double DifferentialCrossSection(const ParticleDef&, const Medium&, double, double);
    double FunctionToDEdxIntegral(const ParticleDef&, const Medium&, double, double);
};

class IonizBergerSeltzerBhabha : public Ionization {
public:
    IonizBergerSeltzerBhabha(const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(const ParticleDef&, const Medium&, double);
    double DifferentialCrossSection(const ParticleDef&, const Medium&, double, double);
    double FunctionToDEdxIntegral(const ParticleDef&, const Medium&, double, double);
};

class IonizBergerSeltzerMoller : public Ionization {
public:
    IonizBergerSeltzerMoller(const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(const ParticleDef&, const Medium&, double);
    double DifferentialCrossSection(const ParticleDef&, const Medium&, double, double);
    double FunctionToDEdxIntegral(const ParticleDef&, const Medium&, double, double);
};

} // namespace PROPOSAL
