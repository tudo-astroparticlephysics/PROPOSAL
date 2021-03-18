
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

#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include "PROPOSAL/methods.h"

#define Q2_PHOTO_PARAM_INTEGRAL_DEC(param)                                     \
    struct Photo##param : public PhotoQ2Integral {                             \
        Photo##param(std::shared_ptr<ShadowEffect>);                           \
        using base_param_t = Photonuclear;                                     \
                                                                               \
        std::unique_ptr<Parametrization<Component>> clone() const final;       \
                                                                               \
        double FunctionToQ2Integral(const ParticleDef&, const Component&,      \
            double energy, double v, double Q2) const;                         \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationName<Photo##param> {                     \
        static constexpr auto value = #param;                                  \
    };                                                                         \
                                                                               \
    template <> struct ParametrizationId<Photo##param> {                       \
        static constexpr size_t value = 1000000005;                            \
    };

namespace PROPOSAL {
class Component;
}

namespace PROPOSAL {
namespace crosssection {
    class ShadowEffect {
    protected:
        size_t hash;

    public:
        ShadowEffect() : hash(0) {};
        virtual ~ShadowEffect() = default;

        virtual double CalculateShadowEffect(
            const Component&, double x, double nu)
            = 0;

        size_t GetHash() const noexcept { return hash; }
    };

    class ShadowDuttaRenoSarcevicSeckel : public ShadowEffect {
    public:
        ShadowDuttaRenoSarcevicSeckel() : ShadowEffect() {
            hash_combine(hash, std::string("duttarenosarcevicseckel"));
        };

        double CalculateShadowEffect(const Component&, double x, double nu);
    };

    class ShadowButkevichMikheyev : public ShadowEffect {
    public:
        ShadowButkevichMikheyev() : ShadowEffect() {
            hash_combine(hash, std::string("butkevichmikheyev"));
        };

        double CalculateShadowEffect(const Component&, double x, double nu);
    };

    class PhotoQ2Integral : public Photonuclear {
    public:
        PhotoQ2Integral(std::shared_ptr<ShadowEffect>);
        virtual ~PhotoQ2Integral() = default;

        virtual double DifferentialCrossSection(const ParticleDef&,
            const Component&, double energy, double v) const;
        virtual double FunctionToQ2Integral(const ParticleDef&,
            const Component&, double energy, double v, double Q2) const = 0;

        std::shared_ptr<ShadowEffect> shadow_effect_;
    };

    Q2_PHOTO_PARAM_INTEGRAL_DEC(AbramowiczLevinLevyMaor91)
    Q2_PHOTO_PARAM_INTEGRAL_DEC(AbramowiczLevinLevyMaor97)
    Q2_PHOTO_PARAM_INTEGRAL_DEC(ButkevichMikheyev)
    Q2_PHOTO_PARAM_INTEGRAL_DEC(RenoSarcevicSu)
    Q2_PHOTO_PARAM_INTEGRAL_DEC(AbtFT)
    Q2_PHOTO_PARAM_INTEGRAL_DEC(BlockDurandHa)

} // namespace crosssection
} // namespace PROPOSAL

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC
