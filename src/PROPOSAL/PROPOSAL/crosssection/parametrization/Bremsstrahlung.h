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

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#define BREMSSTRAHLUNG_DEF(param)                                              \
    struct Brems##param : public Bremsstrahlung {                              \
        Brems##param(bool lpm = false);                                        \
        Brems##param(bool lpm, const ParticleDef&, const Medium&,              \
            double density_correction = 1.0);                                  \
                                                                               \
        std::unique_ptr<Parametrization<Component>> clone() const final;       \
                                                                               \
        double CalculateParametrization(const ParticleDef&, const Component&,  \
            double energy, double v) const final;                              \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationName<Brems##param> {                     \
        static constexpr auto value = #param;                                  \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationId<Brems##param> {                       \
        static constexpr size_t value = 1000000002;                            \
    };

namespace PROPOSAL {
class Interpolant;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace crosssection {
    class BremsLPM;

    class Bremsstrahlung : public Parametrization<Component> {

    protected:
        bool lorenz_;       // enable lorenz cut
        double lorenz_cut_; // in [MeV]
        std::shared_ptr<BremsLPM> lpm_;

    public:
        Bremsstrahlung();

        virtual ~Bremsstrahlung() = default;

        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double, double) const override;
        virtual double CalculateParametrization(
            const ParticleDef&, const Component&, double, double) const = 0;

        double GetLowerEnergyLim(const ParticleDef&) const noexcept final;
        KinematicLimits GetKinematicLimits(
            const ParticleDef&, const Component&, double) const final;
    };

    template <> struct ParametrizationName<Bremsstrahlung> {
        static constexpr auto value = "bremsstrahlung";
    };

    BREMSSTRAHLUNG_DEF(PetrukhinShestakov)
    BREMSSTRAHLUNG_DEF(KelnerKokoulinPetrukhin)
    BREMSSTRAHLUNG_DEF(CompleteScreening)
    BREMSSTRAHLUNG_DEF(AndreevBezrukovBugaev)
    BREMSSTRAHLUNG_DEF(SandrockSoedingreksoRhode)

    class BremsElectronScreening : public Bremsstrahlung {
        std::shared_ptr<Interpolant> interpolant_;

    public:
        BremsElectronScreening(bool lpm = false);
        BremsElectronScreening(bool lpm, const ParticleDef&, const Medium&,
            double density_correction = 1.0);

        std::unique_ptr<Parametrization<Component>> clone() const final;

        double CalculateParametrization(const ParticleDef&, const Component&,
            double energy, double v) const final;
        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double energy, double v) const final;
    };

    template <> struct ParametrizationName<BremsElectronScreening> {
        static constexpr auto value = "ElectronScreening";
    };

    template <> struct ParametrizationId<BremsElectronScreening> {
        static constexpr size_t value = 1000000002;
    };

    // LPM effect object
    class BremsLPM {
        size_t hash;
        double mass_;
        double mol_density_;
        double mass_density_;
        double sum_charge_;
        double density_correction_;
        double eLpm_;

    public:
        BremsLPM(const ParticleDef&, const Medium&, const Bremsstrahlung&,
            double density_correction = 1.0);
        double suppression_factor(
            double energy, double v, const Component&) const;
        size_t GetHash() const noexcept { return hash; }
    };

} // namespace crosssection
} // namespace PROPOSAL

#undef BREMSSTRAHLUNG_DEF
