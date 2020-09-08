
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
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {
namespace crosssection {
    class Ionization : public Parametrization {
    protected:
        EnergyCutSettings cuts_;

        double Delta(const Medium&, double beta, double gamma) const;

    public:
        Ionization(const EnergyCutSettings&);
        virtual ~Ionization() = default;

        using only_stochastic = std::false_type;
        using component_wise = std::false_type;

        virtual double FunctionToDEdxIntegral(
            const ParticleDef&, const Medium&, double, double) const = 0;
        virtual tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Medium&, double) const noexcept
            = 0;
        double FunctionToDNdxIntegral(
            const ParticleDef&, const Medium&, double, double) const;
        double FunctionToDE2dxIntegral(
            const ParticleDef&, const Medium&, double, double) const;
        virtual double DifferentialCrossSection(
            const ParticleDef&, const Medium&, double, double) const = 0;
        double DifferentialCrossSection(
                const ParticleDef&, const Component&, double, double) const final;
        tuple<double, double> GetKinematicLimits(
                const ParticleDef&, const Component&, double) const final;
        double GetLowerEnergyLim(const ParticleDef&) const noexcept override;

        double DifferentialCrossSection(
            const ParticleDef&, const Component&, double, double) const final
        {
            throw std::logic_error("Ionization is not defined component wise.");
        };
        tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Component&, double) const noexcept final
        {
            return std::make_tuple(0.f,0.f);
        };
    };

    class IonizBetheBlochRossi : public Ionization {
        double InelCorrection(
            const ParticleDef&, const Medium&, double, double) const;
        double CrossSectionWithoutInelasticCorrection(
            const ParticleDef&, const Medium&, double, double) const;

    public:
        IonizBetheBlochRossi(const EnergyCutSettings&);
        using base_param_t = Ionization;

        tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Medium&, double) const noexcept final;
        double DifferentialCrossSection(
            const ParticleDef&, const Medium&, double, double) const final;
        double FunctionToDEdxIntegral(
            const ParticleDef&, const Medium&, double, double) const final;
    };

    struct IonizBergerSeltzerBhabha : public Ionization {
        IonizBergerSeltzerBhabha(const EnergyCutSettings&);
        using base_param_t = Ionization;

        tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Medium&, double) const noexcept final;
        double DifferentialCrossSection(
            const ParticleDef&, const Medium&, double, double) const final;
        double FunctionToDEdxIntegral(
            const ParticleDef&, const Medium&, double, double) const final;
    };

    struct IonizBergerSeltzerMoller : public Ionization {
        IonizBergerSeltzerMoller(const EnergyCutSettings&);
        using base_param_t = Ionization;

        tuple<double, double> GetKinematicLimits(
            const ParticleDef&, const Medium&, double) const noexcept final;
        double DifferentialCrossSection(
            const ParticleDef&, const Medium&, double, double) const final;
        double FunctionToDEdxIntegral(
            const ParticleDef&, const Medium&, double, double) const final;
    };
} // crosssection
} // namespace PROPOSAL
