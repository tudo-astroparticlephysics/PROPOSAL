
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
#include "PROPOSAL/crosssection/parametrization/ParametrizationDirect.h"

namespace PROPOSAL {
namespace crosssection {
    class Annihilation : public ParametrizationDirect {
    public:
        Annihilation() = default;

        // no continuous losses
        double CalculatedEdx(double, const ParticleDef&, const Medium&, cut_ptr) override { return 0.; };
        double CalculatedE2dx(double, const ParticleDef&, const Medium&, cut_ptr) override { return 0.; };

        double CalculateCumulativeCrosssection(double, size_t, double, const ParticleDef&, const Medium&, cut_ptr) override {
            throw std::logic_error("No cumulative crosssection defined for Annihilation.");
        };

        double CalculateStochasticLoss(size_t, double, double, const ParticleDef&, const Medium&, cut_ptr) override {
            return 1.; // all energy is always lost, i.e. v=1
        };

        double GetLowerEnergyLim(const ParticleDef&, const Medium&, cut_ptr) const override;

        size_t GetHash(const ParticleDef&, const Medium& m, cut_ptr) const noexcept override;

        InteractionType GetInteractionType() const noexcept override;
    };

    class AnnihilationHeitler : public Annihilation {
    public:
        AnnihilationHeitler();
        std::unique_ptr<ParametrizationDirect> clone() const final;

        double CalculatedNdx(double, const ParticleDef&, const Medium&, cut_ptr) override;
        double CalculatedNdx(double, size_t, const ParticleDef&, const Medium&, cut_ptr) override;
        std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
                double, const ParticleDef&, const Medium&, cut_ptr) override;

    };

    template <> struct ParametrizationName<AnnihilationHeitler> {
        static constexpr auto value = "Heitler";
    };

    template <> struct ParametrizationId<AnnihilationHeitler> {
        static constexpr size_t value = 1000000012;
    };

} // namespace crosssection
} // namespace PROPOSAL
