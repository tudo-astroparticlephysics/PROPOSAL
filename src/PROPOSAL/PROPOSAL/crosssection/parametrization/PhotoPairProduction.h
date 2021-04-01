
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
#include "PROPOSAL/methods.h"

namespace PROPOSAL {
namespace crosssection {
    class PhotoPairProduction : public Parametrization<Component> {
    public:
        PhotoPairProduction() = default;
        virtual ~PhotoPairProduction() = default;

        double GetLowerEnergyLim(const ParticleDef&) const noexcept final;
        KinematicLimits GetKinematicLimits(const ParticleDef&, const Component&,
            double) const noexcept override;
    };

    template <> struct ParametrizationName<PhotoPairProduction> {
        static constexpr auto value = "photopairproduction";
    };

    struct PhotoPairTsai : public PhotoPairProduction {
        PhotoPairTsai() { hash_combine(hash, std::string("tsai")); }
        using base_param_t = PhotoPairProduction;
        std::unique_ptr<Parametrization<Component>> clone() const final;

        virtual double DifferentialCrossSection(
            const ParticleDef&, const Component&, double, double) const;
    };

    template <> struct ParametrizationName<PhotoPairTsai> {
        static constexpr auto value = "Tsai";
    };

    template <> struct ParametrizationId<PhotoPairTsai> {
        static constexpr size_t value = 1000000013;
    };
} // namespace crosssection
} // namespace PROPOSAL
