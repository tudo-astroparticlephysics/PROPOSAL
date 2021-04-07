
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

namespace PROPOSAL {
namespace crosssection {

    template <typename T> struct IonizationDEdxIntegration {
        constexpr static bool value = false;
    };

    struct Ionization : public Parametrization<Medium> {
    protected:
        EnergyCutSettings cuts_;

        double Delta(const Medium&, double beta, double gamma) const;

    public:
        Ionization(const EnergyCutSettings&);
        double GetLowerEnergyLim(const ParticleDef&) const noexcept final;
        virtual ~Ionization() = default;
    };

    template <> struct ParametrizationName<Ionization> {
        static constexpr auto value = "ionization";
    };

    struct IonizBetheBlochRossi : public Ionization {
        double InelCorrection(
            ParticleDef const&, Medium const&, double, double) const;
        double CrossSectionWithoutInelasticCorrection(
            ParticleDef const&, Medium const&, double, double) const;

    public:
        IonizBetheBlochRossi(const EnergyCutSettings&);

        std::unique_ptr<Parametrization<Medium>> clone() const final;

        KinematicLimits GetKinematicLimits(
            ParticleDef const&, Medium const&, double) const final;
        double DifferentialCrossSection(
            ParticleDef const&, Medium const&, double, double) const final;
        double FunctionToDEdxIntegral(
            ParticleDef const&, Medium const&, double, double) const final;
        double IonizationLoss(
                ParticleDef const&, Medium const&, double) const;
    };

    template <> struct ParametrizationName<IonizBetheBlochRossi> {
        static constexpr auto value = "BetheBlochRossi";
    };

    template <> struct ParametrizationId<IonizBetheBlochRossi> {
        static constexpr size_t value = 1000000003;
    };

    template <> struct IonizationDEdxIntegration<IonizBetheBlochRossi> {
        constexpr static bool value = true;
    };

    struct IonizBergerSeltzerBhabha : public Ionization {
        IonizBergerSeltzerBhabha(const EnergyCutSettings&);

        std::unique_ptr<Parametrization<Medium>> clone() const final;

        KinematicLimits GetKinematicLimits(
            const ParticleDef&, const Medium&, double) const final;
        double DifferentialCrossSection(
            const ParticleDef&, const Medium&, double, double) const final;
        double FunctionToDEdxIntegral(
            const ParticleDef&, const Medium&, double, double) const final;
    };

    template <> struct ParametrizationName<IonizBergerSeltzerBhabha> {
        static constexpr auto value = "BergerSeltzerBhabha";
    };

    template <> struct ParametrizationId<IonizBergerSeltzerBhabha> {
        static constexpr size_t value = 1000000003;
    };

    struct IonizBergerSeltzerMoller : public Ionization {
        IonizBergerSeltzerMoller(const EnergyCutSettings&);

        std::unique_ptr<Parametrization<Medium>> clone() const final;

        KinematicLimits GetKinematicLimits(
            const ParticleDef&, const Medium&, double) const final;
        double DifferentialCrossSection(
            const ParticleDef&, const Medium&, double, double) const final;
        double FunctionToDEdxIntegral(
            const ParticleDef&, const Medium&, double, double) const final;
    };

    template <> struct ParametrizationName<IonizBergerSeltzerMoller> {
        static constexpr auto value = "BergerSeltzerMoller";
    };

    template <> struct ParametrizationId<IonizBergerSeltzerMoller> {
        static constexpr size_t value = 1000000003;
    };

} // crosssection
} // namespace PROPOSAL
