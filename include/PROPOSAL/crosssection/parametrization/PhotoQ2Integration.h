
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

#define Q2_PHOTO_PARAM_INTEGRAL_DEC(param)                                     \
    struct Photo##param : public PhotoQ2Integral {                             \
        Photo##param(std::shared_ptr<ShadowEffect>);                           \
        using base_param_t = Photonuclear;                                     \
        double FunctionToQ2Integral(const ParticleDef&, const Component&,      \
            double energy, double v, double Q2) const;                         \
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
        ShadowEffect() = default;
        virtual ~ShadowEffect() = default;

        virtual double CalculateShadowEffect(
            const Component&, double x, double nu)
            = 0;

        size_t GetHash() const noexcept { return hash; }
    };

    class ShadowDuttaRenoSarcevicSeckel : public ShadowEffect {
    public:
        ShadowDuttaRenoSarcevicSeckel() = default;

        double CalculateShadowEffect(const Component&, double x, double nu);
    };

    class ShadowButkevichMikheyev : public ShadowEffect {
    public:
        ShadowButkevichMikheyev() = default;

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

    /* // Factory pattern functions */

    /* template <typename P, typename M> */
    /* using photoQ2_func_ptr */
    /*     = cross_t_ptr<P, M> (*)(P, M, std::shared_ptr<const
     * EnergyCutSettings>, */
    /*         std::shared_ptr<ShadowEffect>, bool); */

    /* template <typename Param, typename P, typename M> */
    /* cross_t_ptr<P, M> create_photoQ2(P p_def, M medium, */
    /*     std::shared_ptr<const EnergyCutSettings> cuts, */
    /*     std::shared_ptr<ShadowEffect> shadow, bool interpol) */
    /* { */
    /*     auto param = Param(shadow); */
    /*     return make_crosssection(param, p_def, medium, cuts, interpol); */
    /* } */

    /* template <typename P, typename M> */
    /* static std::map<std::string, photoQ2_func_ptr<P, M>> photoQ2_map = { */
    /*     { "abramowiczlevinlevymaor91", */
    /*         create_photoQ2<PhotoAbramowiczLevinLevyMaor91, P, M> }, */
    /*     { "abramowiczlevinlevymaor97", */
    /*         create_photoQ2<PhotoAbramowiczLevinLevyMaor97, P, M> }, */
    /*     { "butkevichmikheyev", create_photoQ2<PhotoButkevichMikheyev, P, M>
     * }, */
    /*     { "renosarcevicsu", create_photoQ2<PhotoRenoSarcevicSu, P, M> } */
    /* }; */

    /* using shadow_func_ptr = std::shared_ptr<ShadowEffect> (*)(); */

    /* template <typename Param> std::shared_ptr<ShadowEffect> create_shadow()
     */
    /* { */
    /*     return std::make_shared<Param>(); */
    /* } */

    /* static std::map<std::string, shadow_func_ptr> shadow_map = { */
    /*     { "duttarenosarcevicseckel", */
    /*         create_shadow<ShadowDuttaRenoSarcevicSeckel> }, */
    /*     { "butkevichmikheyev", create_shadow<ShadowButkevichMikheyev> }, */
    /* }; */

    /* template <typename P, typename M> */
    /* cross_t_ptr<P, M> make_photonuclearQ2(P p_def, M medium, */
    /*     std::shared_ptr<const EnergyCutSettings> cuts, bool interpol, */
    /*     const std::string& param_name, const std::string& shadow_name) */
    /* { */
    /*     std::string name = param_name; */
    /*     std::transform( */
    /*         param_name.begin(), param_name.end(), name.begin(), ::tolower);
     */
    /*     auto it = photoQ2_map<P, M>.find(name); */
    /*     if (it == photoQ2_map<P, M>.end()) */
    /*         throw std::invalid_argument( */
    /*             "Unknown parametrization for photonuclear"); */

    /*     std::string name2 = shadow_name; */
    /*     std::transform( */
    /*         shadow_name.begin(), shadow_name.end(), name2.begin(),
     * ::tolower); */
    /*     auto it_shadow = shadow_map.find(name2); */
    /*     if (it_shadow == shadow_map.end()) */
    /*         throw std::logic_error("Shadow effect name unknown"); */
    /*     return it->second(p_def, medium, cuts, it_shadow->second(),
     * interpol); */
    /* } */

    /* template <typename P, typename M> */
    /* cross_t_ptr<P, M> make_photonuclearQ2(P p_def, M medium, */
    /*     std::shared_ptr<const EnergyCutSettings> cuts, bool interpol, */
    /*     const nlohmann::json& config) */
    /* { */
    /*     if (!config.contains("parametrization")) */
    /*         throw std::logic_error( */
    /*             "No parametrization passed for photonuclear"); */

    /*     std::string param_name = config["parametrization"]; */
    /*     std::string shadow_name = config.value("shadow",
     * "ButkevichMikheyev"); */

    /*     return make_photonuclearQ2( */
    /*         p_def, medium, cuts, interpol, param_name, shadow_name); */
    /* } */

} // namespace crosssection
} // namespace PROPOSAL

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC
