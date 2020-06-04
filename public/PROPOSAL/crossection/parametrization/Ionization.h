
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
protected:
    EnergyCutSettings cuts_;

    double Delta(const Medium&, double beta, double gamma);

public:
    Ionization(const EnergyCutSettings&);
    virtual ~Ionization() = default;

    using only_stochastic = std::false_type;
    using component_wise = std::false_type;

    virtual double FunctionToDEdxIntegral(
        const ParticleDef&, const Medium&, double, double)
        = 0;
    virtual tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Medium&, double) const noexcept
        = 0;
    double FunctionToDNdxIntegral(
        const ParticleDef&, const Medium&, double, double);
    double FunctionToDE2dxIntegral(
        const ParticleDef&, const Medium&, double, double);
    virtual double DifferentialCrossSection(
        const ParticleDef&, const Medium&, double, double)
        = 0;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
};

class IonizBetheBlochRossi : public Ionization {
    double InelCorrection(const ParticleDef&, const Medium&, double, double);
    double CrossSectionWithoutInelasticCorrection(
        const ParticleDef&, const Medium&, double, double);

public:
    IonizBetheBlochRossi(const EnergyCutSettings&);
    using base_param_t = Ionization;

    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Medium&, double) const noexcept override;
    double DifferentialCrossSection(
        const ParticleDef&, const Medium&, double, double) override;
    double FunctionToDEdxIntegral(
        const ParticleDef&, const Medium&, double, double) override;
};

struct IonizBergerSeltzerBhabha : public Ionization {
    IonizBergerSeltzerBhabha(const EnergyCutSettings&);
    using base_param_t = Ionization;

    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Medium&, double) const noexcept override;
    double DifferentialCrossSection(
        const ParticleDef&, const Medium&, double, double) override;
    double FunctionToDEdxIntegral(
        const ParticleDef&, const Medium&, double, double) override;
};

struct IonizBergerSeltzerMoller : public Ionization {
    IonizBergerSeltzerMoller(const EnergyCutSettings&);
    using base_param_t = Ionization;

    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Medium&, double) const noexcept override;
    double DifferentialCrossSection(
        const ParticleDef&, const Medium&, double, double) override;
    double FunctionToDEdxIntegral(
        const ParticleDef&, const Medium&, double, double) override;
};

} // namespace PROPOSAL
