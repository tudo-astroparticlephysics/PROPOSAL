
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
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"

namespace PROPOSAL {
class Ionization : public Parametrization {
protected:
    EnergyCutSettings cuts_;
    double X0_;
    double X1_;
    double D0_;
    double m_;
    double C_;
    double a_;
    double I_;
    double ZA_;

    double Delta(double beta, double gamma);

public:
    Ionization(const ParticleDef&, const Medium&, const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(double energy) = 0;
    double DifferentialCrossSection(double energy, double v) = 0;
    double FunctionToDEdxIntegral(double energy, double v) = 0;
};

class IonizBetheBlochRossi : public Ionization {
    double InelCorrection(double energy, double v);
    double CrossSectionWithoutInelasticCorrection(double energy, double v);

public:
    IonizBetheBlochRossi(
        const ParticleDef&, const Medium&, const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(double energy) override;
    double DifferentialCrossSection(double energy, double v) override;
    double FunctionToDEdxIntegral(double energy, double v) override;
};

class IonizBergerSeltzerBhabha : public Ionization {
public:
    IonizBergerSeltzerBhabha(
        const ParticleDef&, const Medium&, const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(double energy) override;
    double DifferentialCrossSection(double energy, double v) override;
    double FunctionToDEdxIntegral(double energy, double v) override;
};

class IonizBergerSeltzerMoller : public Ionization {
public:
    IonizBergerSeltzerMoller(
        const ParticleDef&, const Medium&, const EnergyCutSettings&);

    KinematicLimits GetKinematicLimits(double energy) override;
    double DifferentialCrossSection(double energy, double v) override;
    double FunctionToDEdxIntegral(double energy, double v) override;
};

} // namespace PROPOSAL
