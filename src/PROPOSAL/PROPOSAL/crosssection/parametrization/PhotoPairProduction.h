
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

#include <string>
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {
namespace crosssection {
    class PhotoPairLPM;

    class PhotoPairProduction : public Parametrization<Component> {
    protected:
        std::shared_ptr<PhotoPairLPM> lpm_;
        double density_correction_; // correction to standard medium density for LPM
    public:
        PhotoPairProduction();
        virtual ~PhotoPairProduction() = default;

        double GetLowerEnergyLim(const ParticleDef&) const noexcept final;
        KinematicLimits GetKinematicLimits(const ParticleDef&, const Component&,
            double) const noexcept override;
    };

    template <> struct ParametrizationName<PhotoPairProduction> {
        static constexpr auto value = "photopairproduction";
    };

    struct PhotoPairTsai : public PhotoPairProduction {
        PhotoPairTsai(bool lpm = false);
        PhotoPairTsai(bool lpm, const ParticleDef&, const Medium&,
                          double density_correction = 1.0);
        using base_param_t = PhotoPairProduction;
        std::unique_ptr<Parametrization<Component>> clone() const final;

        virtual double DifferentialCrossSection(
            const ParticleDef&, const Component&, double, double) const;
    };

    struct PhotoPairKochMotz : public PhotoPairProduction {
        PhotoPairKochMotz(bool lpm = false);
        PhotoPairKochMotz(bool lpm, const ParticleDef&, const Medium&,
                          double density_correction = 1.0);
        using base_param_t = PhotoPairProduction;
        std::unique_ptr<Parametrization<Component>> clone() const final;

        virtual double DifferentialCrossSection(
                const ParticleDef&, const Component&, double, double) const;

        private:
            std::shared_ptr<Interpolant> interpolant_;
            double DifferentialCrossSectionWithoutA(
                    const ParticleDef&, const Component&, double, double) const;
    };

    template <> struct ParametrizationName<PhotoPairTsai> {
        static constexpr auto value = "Tsai";
    };

    template <> struct ParametrizationName<PhotoPairKochMotz> {
        static constexpr auto value = "kochmotz";
    };

    template <> struct ParametrizationId<PhotoPairTsai> {
        static constexpr size_t value = 1000000013;
    };

    template <> struct ParametrizationId<PhotoPairKochMotz> {
        static constexpr size_t value = 1000000013;
    };

    // LPM effect object
    class PhotoPairLPM {
        size_t hash;
        double mol_density_;
        double mass_density_;
        double sum_charge_;
        double eLpm_;

    public:
        PhotoPairLPM(const ParticleDef&, const Medium&, const PhotoPairProduction&);
        double suppression_factor(double energy, double x, const Component&,
                                  double density_correction = 1.0) const;
        size_t GetHash() const noexcept { return hash; }
    };

} // namespace crosssection
} // namespace PROPOSAL
